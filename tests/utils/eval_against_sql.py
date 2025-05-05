#!/usr/bin/env python3

"""Module for evaluating fresh data against cached data for testing."""

import os
from contextlib import contextmanager
from typing import Dict, List

import numpy as np
import pandas as pd
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
) -> Dict[int, pd.DataFrame]:
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
    Dict[int, pd.DataFrame]
        Dictionary mapping shot IDs to retrieved fresh data.
    """
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        time_setting="disruption_warning",
        run_columns=test_columns,
        only_requested_columns=True,
    )
    return get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dict",
        log_settings=LogSettings(
            file_path=log_file_path,
            file_level="DEBUG",
            console_level=console_log_level,
        ),
    )


def get_cached_from_fresh(
    tokamak: Tokamak,
    shotlist: List[int],
    fresh_data: Dict[int, pd.DataFrame],
    test_columns: List[str] = None,
) -> Dict[int, pd.DataFrame]:
    """
    Get cached SQL data for a list of shots and map onto the timebase of the supplied
    fresh MDSplus data.

    Returns
    -------
    Dict[int, pd.DataFrame]
        Dictionary mapping shot IDs to retrieved SQL data.
    """
    # Mapping SQL data onto the MDSplus timebase means SQL data needs time data
    merge_col = "time"
    if test_columns is None:
        test_columns = ["*"]
    elif merge_col not in test_columns:
        test_columns.append(merge_col)

    db = ShotDatabase.from_config(tokamak=tokamak)
    shot_data = {}
    for shot_id in shotlist:
        times = fresh_data[shot_id]["time"]
        sql_data = db.get_shots_data([shot_id], cols=test_columns)

        if sql_data.empty:
            shot_data[shot_id] = pd.DataFrame()
            continue

        shot_data[shot_id] = pd.merge_asof(
            times.to_frame(),
            sql_data,
            on=merge_col,
            direction="nearest",
            tolerance=config().time.time_const,
        )
    return shot_data


def eval_shots_against_cache(
    shotlist: List[int],
    fresh_data: Dict[int, pd.DataFrame],
    cache_data: Dict[int, pd.DataFrame],
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
            fresh_shot_data, cache_shot_data = fresh_data[shot_id], cache_data[shot_id]
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
    fresh_shot_data: pd.DataFrame,
    cache_shot_data: pd.DataFrame,
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
        fresh_column_data=fresh_shot_data.get(data_column, None),
        cache_column_data=cache_shot_data.get(data_column, None),
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
            "column {data_difference.data_column}"
        )

    return data_difference


def eval_against_cache(
    tokamak: Tokamak,
    shotlist: List[int],
    expected_failure_columns: List[str],
    test_columns=None,
    console_log_level="WARNING",
) -> Dict[int, pd.DataFrame]:
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
    Dict[int, pd.DataFrame]
        A dictionary mapping shot identifiers to their corresponding data
        differences as DataFrames.
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
