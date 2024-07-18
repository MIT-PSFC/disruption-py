#!/usr/bin/env python3

import inspect
import logging
import os
import time
from contextlib import contextmanager
from tempfile import mkdtemp
from typing import Callable, Dict, List

import numpy as np
import pandas as pd

from disruption_py.config import config
from disruption_py.core.utils.math import matlab_gradient_1d_vectorized
from disruption_py.io.sql import ShotDatabase
from disruption_py.machine import tokamak
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.data_difference import DataDifference


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
        run_tags=[] if test_columns else None,
        run_columns=test_columns if test_columns else [],
        only_requested_columns=test_columns,
    )
    shot_data = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting="dict",
        log_settings=LogSettings(
            log_to_console=False,
            log_file_path=log_file_path,
            log_file_write_mode="w",
            file_log_level=logging.DEBUG,
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
    fail_quick: bool = False,
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
                fail_quick=fail_quick,
                expect_failure=expect_failure,
            )
            data_differences.append(data_difference)
    return data_differences


def eval_shot_against_sql(
    shot_id: int,
    mdsplus_shot_data: pd.DataFrame,
    sql_shot_data: pd.DataFrame,
    data_column: str,
    fail_quick: bool = False,
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
        missing_mdsplus_data=missing_mdsplus_data,
        missing_sql_data=missing_sql_data,
        expect_failure=expect_failure,
    )
    if fail_quick and not (missing_mdsplus_data or missing_sql_data):
        expectation = "failure" if expect_failure else "success"
        failure = "failed" if data_difference.failed else "succeeded"
        assert (
            not data_difference.failed
        ), "Expected {expectation} and {failure}:\n{report}".format(
            expectation=expectation,
            failure=failure,
            report=data_difference.column_mismatch_string,
        )

    return data_difference


def get_failure_statistics_string(
    data_differences: list["DataDifference"], data_column=None
):
    data_difference_by_column = {}
    for data_difference in data_differences:
        data_difference_by_column.setdefault(data_difference.data_column, []).append(
            data_difference
        )

    failure_strings = {}
    failed_columns, succeeded_columns, missing_data_columns = set(), set(), set()
    matches_expected_failures_columns, not_matches_expected_failures_columns = (
        set(),
        set(),
    )
    for ratio_data_column, column_data_differences in data_difference_by_column.items():
        failures = [
            data_difference.shot_id
            for data_difference in column_data_differences
            if data_difference.failed
        ]
        failed = len(failures) > 0

        all_missing_data = all(
            [
                data_difference.missing_data
                for data_difference in column_data_differences
            ]
        )

        anomaly_count = sum(
            [
                data_difference.num_anomalies
                for data_difference in column_data_differences
            ]
        )
        timebase_count = sum(
            [
                data_difference.timebase_length
                for data_difference in column_data_differences
            ]
        )

        expect_failure = any(
            [
                data_difference.expect_failure
                for data_difference in column_data_differences
            ]
        )
        matches_expected_failure = expect_failure == failed

        # failure string
        failure_string_lines = [
            f"Column {ratio_data_column} {'FAILED' if failed else 'SUCCEEDED'}",
            f"Matches expected failures: {matches_expected_failure}",
            f"Total Entry Failure Rate: {anomaly_count / timebase_count * 100:.2f}%",
        ]
        failure_string = "\n".join(failure_string_lines)

        # condition string
        conditions: Dict[str, Callable[[DataDifference], bool]] = {
            "Shots expected to fail that failed": lambda data_difference: data_difference.expect_failure
            and data_difference.failed,
            "Shots expected to succeed that failed": lambda data_difference: not data_difference.expect_failure
            and data_difference.failed,
            "Shots expected to fail that succeeded": lambda data_difference: data_difference.expect_failure
            and not data_difference.failed,
            "Shots expected to succeed that succeeded": lambda data_difference: not data_difference.expect_failure
            and not data_difference.failed,
            "Shots missing sql data": lambda data_difference: data_difference.missing_sql_data,
            "Shots missing mdsplus data": lambda data_difference: data_difference.missing_mdsplus_data,
        }
        condition_results = {}
        for condition_name, condition in conditions.items():
            shotlist = [
                data_difference.shot_id
                for data_difference in column_data_differences
                if condition(data_difference)
            ]
            if len(shotlist) > 0:
                condition_results[condition_name] = shotlist
        condition_string = "\n".join(
            [
                f"{condition_name} ({len(condition_result)} shots): {condition_result}"
                for condition_name, condition_result in condition_results.items()
            ]
        )

        # combine the string parts together
        failure_strings[ratio_data_column] = failure_string + "\n" + condition_string

        if all_missing_data:
            missing_data_columns.add(ratio_data_column)
        elif failed:
            failed_columns.add(ratio_data_column)
        else:
            succeeded_columns.add(ratio_data_column)

        if matches_expected_failure:
            matches_expected_failures_columns.add(ratio_data_column)
        else:
            not_matches_expected_failures_columns.add(ratio_data_column)

    if data_column is not None:
        return failure_strings.get(data_column, "")
    else:
        summary_string = f"""\
        ___________________________________________________________________________________________________
        SUMMARY
        Columns with a failure:
        {"None" if len(failed_columns) == 0 else ""}{", ".join(failed_columns)}
        
        Columns without a failure:
        {"None" if len(succeeded_columns) == 0 else ""}{", ".join(succeeded_columns)}
        
        Columns lacking data for comparison from sql or mdsplus sources:
        {"None" if len(missing_data_columns) == 0 else ""}{", ".join(missing_data_columns)}
        ___________________________________________________________________________________________________
        
        Columns that match expected failures:
        {"None" if len(matches_expected_failures_columns) == 0 else ""}{", ".join(matches_expected_failures_columns)}
        
        Columns that do not match expected failures:
        {"None" if len(not_matches_expected_failures_columns) == 0 else ""}{", ".join(not_matches_expected_failures_columns)}
        """
        return (
            "\n\n".join(failure_strings.values())
            + "\n\n"
            + inspect.cleandoc(summary_string)
        )


def eval_against_sql(
    tokamak: Tokamak,
    shotlist: List[int],
    expected_failure_columns: List[str],
    fail_quick: bool,
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
            tokamak,
            shotlist,
            os.path.join(tempfolder, "data_retrieval.log"),
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
        fail_quick=fail_quick,
        expected_failure_columns=expected_failure_columns,
    )

    return data_differences
