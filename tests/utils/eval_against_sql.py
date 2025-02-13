#!/usr/bin/env python3

"""Module for evaluating fresh data against cached data for testing."""

import os
from contextlib import contextmanager
from typing import Dict, List

import numpy as np
import xarray as xr
from loguru import logger

from disruption_py.config import config
from disruption_py.core.utils.math import matlab_gradient_1d_vectorized
from disruption_py.core.utils.misc import get_temporary_folder
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.data_difference import DataDifference
from tests.utils.factory import get_tokamak_test_columns


def get_fresh_data(
    tokamak: Tokamak,
    shotlist: List[int],
    log_file_path: str,
    test_columns: List[str] = None,
    console_log_level: str = "WARNING",
) -> xr.Dataset:
    """
    Get fresh data for a list of shots.

    Parameters
    ----------
    tokamak : Tokamak
        The tokamak for which to retrieve data.
    shotlist : List[int]
        A list of shot IDs to retrieve data for.
    log_file_path : str
        The path to the log file.
    test_columns : List[str], optional
        A list of columns to retrieve.
    console_log_level : str, optional
        The log level for console output. Default is "WARNING".

    Returns
    -------
    xr.Dataset
        Fresh dataset.
    """
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        time_setting="disruption_warning",
        run_columns=test_columns,
        only_requested_columns=True,
    )
    shot_data = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dataset",
        log_settings=LogSettings(
            log_to_console=True,
            log_file_path=log_file_path,
            log_file_write_mode="w",
            file_log_level="DEBUG",
            console_log_level=console_log_level,
        ),
    )
    return shot_data


def get_cached_from_fresh(
    tokamak: Tokamak,
    shotlist: List[int],
    fresh_data: xr.Dataset,
    test_columns: List[str] = None,
) -> xr.Dataset:
    """
    Get cached SQL data for a list of shots and map onto the timebase of the supplied
    fresh MDSplus data.

    Returns
    -------
    xr.Dataset
        Dataset of retrieved SQL data.
    """
    # Mapping SQL data onto the MDSplus timebase means SQL data needs time data
    merge_col = "time"
    if test_columns is None:
        test_columns = ["*"]
    elif merge_col not in test_columns:
        test_columns.append(merge_col)

    index = ["shot", "time"]

    db = ShotDatabase.from_config(tokamak=tokamak)
    shot_data_list = []
    for shot_id in shotlist:
        sql_data_df = db.get_shots_data([shot_id], cols=test_columns).astype(float)
        if sql_data_df.empty:
            shot_data_list.append(xr.Dataset(coords={"shot": [shot_id], "time": []}))
            continue
        sql_data_df["shot"] = shot_id

        # Some shots, like 1150805012 on C-Mod, have multiple entries for the same
        # shot and time. You can't create a Dataset from a DF with duplicate indices
        sql_data_df.drop_duplicates(subset=index, inplace=True)
        sql_data = xr.Dataset.from_dataframe(sql_data_df.set_index(index))
        shot_data_list.append(sql_data)
    shot_data = xr.concat(shot_data_list, dim="shot")
    shot_data = shot_data.reindex(
        {"time": fresh_data.time}, method="nearest", tolerance=config().time_const
    )
    return shot_data


def eval_shots_against_cache(
    shotlist: List[int],
    fresh_data: xr.Dataset,
    cache_data: xr.Dataset,
    data_columns: List[str],
    expected_failure_columns: List[str] = None,
) -> List["DataDifference"]:
    """
    Test if the difference between the two data sources is within tolerance.
    """
    if expected_failure_columns is None:
        expected_failure_columns = []

    data_differences: List[DataDifference] = []
    for data_column in data_columns:
        for shot_id in shotlist:
            fresh_shot_data, cache_shot_data = fresh_data.sel(
                shot=shot_id
            ), cache_data.sel(shot=shot_id)
            expect_failure = data_column in expected_failure_columns

            data_difference = eval_shot_against_cache(
                shot_id=shot_id,
                fresh_shot_data=fresh_shot_data,
                cache_shot_data=cache_shot_data,
                data_column=data_column,
                expect_failure=expect_failure,
            )
            data_differences.append(data_difference)
    return data_differences


def eval_shot_against_cache(
    shot_id: int,
    fresh_shot_data: xr.Dataset,
    cache_shot_data: xr.Dataset,
    data_column: str,
    expect_failure: bool = False,
) -> "DataDifference":
    """
    Test if the difference between the two data is within tolerance.
    """
    missing_fresh_data = data_column not in fresh_shot_data
    missing_cache_data = data_column not in cache_shot_data
    data_difference = DataDifference(
        shot_id=shot_id,
        data_column=data_column,
        fresh_data_array=fresh_shot_data.get(data_column, None),
        cache_data_array=cache_shot_data.get(data_column, None),
        fresh_time=fresh_shot_data.get("time"),
        cache_time=cache_shot_data.get("time"),
        missing_fresh_data=missing_fresh_data,
        missing_cache_data=missing_cache_data,
        expect_failure=expect_failure,
    )

    # Condition on both failing and expecting to fail to log expected & unexpected
    # failures and unexpected successes
    if data_difference.failed or data_difference.expect_failure:
        expectation = "failure" if expect_failure else "success"
        failure = "failed" if data_difference.failed else "succeeded"
        logger.debug(
            "Expected {expectation} and {failure}:\n{mismatch_string}",
            expectation=expectation,
            failure=failure,
            mismatch_string=data_difference.column_mismatch_string,
        )
    # Python tests should not assert for expected failures to only catch unexpected failures
    # Pytest should assert for expected failures to confirm the test fails or
    # to catch unexpected successes
    if "PYTEST_CURRENT_TEST" in os.environ or not expect_failure:
        assert not data_difference.failed, (
            f"Comparison failed on shot {data_difference.shot_id}, "
            f"column {data_difference.data_column}"
        )

    return data_difference


def eval_against_cache(
    tokamak: Tokamak,
    shotlist: List[int],
    expected_failure_columns: List[str],
    test_columns=None,
    console_log_level="WARNING",
) -> Dict[int, xr.Dataset]:
    """
    Evaluate fresh data against cached data for specified shots.

    This function retrieves fresh data from a tokamak and compares it against
    cached data, identifying any differences. It temporarily patches the NumPy
    gradient function to ensure compatibility with MATLAB-style gradient calculations.

    Parameters
    ----------
    tokamak : Tokamak
        The tokamak object used to retrieve data.
    shotlist : List[int]
        A list of shot identifiers to evaluate.
    expected_failure_columns : List[str]
        A list of columns that are expected to fail during evaluation.
    test_columns : List[str], optional
        A list of columns to test against the cached data. If None, the function
        will determine the columns based on the available data.
    console_log_level : str, optional
        The log level for console output. Default is "WARNING".

    Returns
    -------
    List[DataDifference]
        A list of DataDifference objects; one for each shot.
    """

    tempfolder = get_temporary_folder()
    print(f"Outputting to temporary folder: {tempfolder}")

    @contextmanager
    def monkey_patch_numpy_gradient():
        original_function = np.gradient
        np.gradient = matlab_gradient_1d_vectorized
        try:
            yield
        finally:
            np.gradient = original_function

    with monkey_patch_numpy_gradient():
        fresh_data = get_fresh_data(
            tokamak=tokamak,
            shotlist=shotlist,
            log_file_path=os.path.join(tempfolder, "data_retrieval.log"),
            test_columns=test_columns,
            console_log_level=console_log_level,
        )
    cache_data = get_cached_from_fresh(tokamak, shotlist, fresh_data, test_columns)

    if test_columns is None:
        test_columns = get_tokamak_test_columns(tokamak)

    data_differences = eval_shots_against_cache(
        shotlist=shotlist,
        fresh_data=fresh_data,
        cache_data=cache_data,
        data_columns=test_columns,
        expected_failure_columns=expected_failure_columns,
    )

    return data_differences
