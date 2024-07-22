#!/usr/bin/env python3

import inspect
import logging
import os
import time
from contextlib import contextmanager
from tempfile import mkdtemp
from typing import Dict, List

import numpy as np
import pandas as pd
import pytest

from disruption_py.config import config
from disruption_py.core.utils.math import matlab_gradient_1d_vectorized
from disruption_py.io.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.data_difference import DataDifference

logger = logging.getLogger("disruption_py")


def get_mdsplus_data(
    tokamak: Tokamak,
    shotlist: List[int],
    log_file_path: str,
    test_columns: List[str] = None,
) -> Dict[int, pd.DataFrame]:
    """
    Get MDSplus data for a list of shots.

     Returns
    -------
    Dict[int, pd.DataFrame]
        Dictionary mapping shot IDs to retrieved MDSplus data.
    """
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        time_setting="disruption_warning",
        run_tags=[] if test_columns else ["all"],
        run_columns=test_columns if test_columns else [],
        only_requested_columns=test_columns,
    )
    shot_data = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dict",
        log_settings=LogSettings(
            log_to_console=True,
            log_file_path=log_file_path,
            log_file_write_mode="w",
            file_log_level=logging.DEBUG,
            console_log_level=logging.DEBUG,
        ),
    )
    return shot_data


def get_sql_data_for_mdsplus(
    tokamak: Tokamak,
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    test_columns: List[str] = None,
) -> Dict[int, pd.DataFrame]:
    """
    Get SQL data for a list of shots and map onto the timebase of the supplied MDSplus data.

    Returns
    -------
    Dict[int, pd.DataFrame]
        Dictionary mapping shot IDs to retrieved SQL data.
    """
    # Mapping SQL data onto the MDSplus timebase means SQL data needs time data
    MERGE_COL = "time"
    if test_columns is None:
        test_columns = ["*"]
    elif MERGE_COL not in test_columns:
        test_columns.append(MERGE_COL)

    db = ShotDatabase.from_config(tokamak=tokamak)
    shot_data = {}
    for shot_id in shotlist:
        times = mdsplus_data[shot_id]["time"]
        sql_data = db.get_shots_data([shot_id], cols=test_columns)
        shot_data[shot_id] = pd.merge_asof(
            times.to_frame(),
            sql_data,
            on=MERGE_COL,
            direction="nearest",
            tolerance=config().TIME_CONST,
        )
    return shot_data


def eval_shots_against_sql(
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    sql_data: Dict[int, pd.DataFrame],
    data_columns: List[str],
    expected_failure_columns: List[str] = None,
) -> List["DataDifference"]:
    """
    Test if the difference between the two data is within tolerance.
    """
    if expected_failure_columns is None:
        expected_failure_columns = []

    data_differences: List[DataDifference] = []
    for data_column in data_columns:
        for shot_id in shotlist:
            mdsplus_shot_data, sql_shot_data = mdsplus_data[shot_id], sql_data[shot_id]
            expect_failure = data_column in expected_failure_columns

            data_difference = eval_shot_against_sql(
                shot_id=shot_id,
                mdsplus_shot_data=mdsplus_shot_data,
                sql_shot_data=sql_shot_data,
                data_column=data_column,
                expect_failure=expect_failure,
            )
            data_differences.append(data_difference)
    return data_differences


def eval_shot_against_sql(
    shot_id: int,
    mdsplus_shot_data: pd.DataFrame,
    sql_shot_data: pd.DataFrame,
    data_column: str,
    expect_failure: bool = False,
) -> "DataDifference":
    """
    Test if the difference between the two data is within tolerance.
    """
    missing_mdsplus_data = data_column not in mdsplus_shot_data
    missing_sql_data = data_column not in sql_shot_data
    data_difference = DataDifference(
        shot_id=shot_id,
        data_column=data_column,
        mdsplus_column_data=mdsplus_shot_data.get(data_column, None),
        sql_column_data=sql_shot_data.get(data_column, None),
        mds_time=mdsplus_shot_data["time"],
        sql_time=sql_shot_data["time"],
        missing_mdsplus_data=missing_mdsplus_data,
        missing_sql_data=missing_sql_data,
        expect_failure=expect_failure,
    )

    # Condition on both failing and expecting to fail to log expected & unexpected
    # failures and unexpected successes
    if data_difference.failed or data_difference.expect_failure:
        expectation = "failure" if expect_failure else "success"
        failure = "failed" if data_difference.failed else "succeeded"
        logger.debug(
            "Expected {expectation} and {failure}:\n{report}".format(
                expectation=expectation,
                failure=failure,
                report=data_difference.column_mismatch_string,
            )
        )
    assert not data_difference.failed, "Comparison failed"

    return data_difference


def eval_against_sql(
    tokamak: Tokamak,
    shotlist: List[int],
    expected_failure_columns: List[str],
    test_columns=None,
) -> Dict[int, pd.DataFrame]:

    tempfolder = mkdtemp(prefix=f"disruptionpy-{time.strftime('%y%m%d-%H%M%S')}-")
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
        mdsplus_data = get_mdsplus_data(
            tokamak=tokamak,
            shotlist=shotlist,
            log_file_path=os.path.join(tempfolder, "data_retrieval.log"),
            test_columns=test_columns,
        )
    sql_data = get_sql_data_for_mdsplus(tokamak, shotlist, mdsplus_data, test_columns)

    if test_columns is None:
        mdsplus_columns = set().union(*(df.columns for df in mdsplus_data.values()))
        sql_columns = set().union(*(df.columns for df in sql_data.values()))
        test_columns = sorted(mdsplus_columns.intersection(sql_columns))

    data_differences = eval_shots_against_sql(
        shotlist=shotlist,
        mdsplus_data=mdsplus_data,
        sql_data=sql_data,
        data_columns=test_columns,
        expected_failure_columns=expected_failure_columns,
    )

    return data_differences
