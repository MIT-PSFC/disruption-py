#!/usr/bin/env python3

"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""

import argparse
import re
from typing import Dict, List

import pandas as pd
import pytest

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from tests.utils.eval_against_sql import (
    eval_against_sql,
    eval_shots_against_sql,
    get_failure_statistics_string,
    get_mdsplus_data,
    get_sql_data_for_mdsplus,
)
from tests.utils.factory import (
    get_tokamak_test_expected_failure_columns,
    get_tokamak_test_shotlist,
)


def extract_param(config):
    """Extract the data column from the pytest command.

    E.g. will return ip given
    `pytest -s tests/test_against_sql.py -k test_data_columns[ip]`

    Params:
        config: pytestconfig fixture

    Returns:
        List[str], the data column if it exists, otherwise None.
    """
    args = config.invocation_params.args
    m = re.search(r"\[(.+)\]$", args[-1])
    param = [m.group(1)] if m is not None else None
    return param


@pytest.fixture(scope="module")
def mdsplus_data(
    tokamak: Tokamak, shotlist: List[int], module_file_path_f, pytestconfig
) -> Dict[int, pd.DataFrame]:
    return get_mdsplus_data(
        tokamak=tokamak,
        shotlist=shotlist,
        log_file_path=module_file_path_f(".log"),
        test_columns=extract_param(pytestconfig),
    )


@pytest.fixture(scope="module")
def sql_data(
    tokamak: Tokamak,
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    return get_sql_data_for_mdsplus(
        tokamak=tokamak,
        shotlist=shotlist,
        mdsplus_data=mdsplus_data,
        test_columns=extract_param(pytestconfig),
    )


def test_data_columns(
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    sql_data: Dict[int, pd.DataFrame],
    data_column,
    expected_failure_columns: List[str],
    fail_quick: bool,
):
    """
    Test that the data columns are the same between MDSplus and SQL across specified data columns.

    Data column is parameterized in pytest_generate_tests.
    """
    # if data_column in expected_failure_columns:
    #     request.node.add_marker(pytest.mark.xfail(reason='column expected failure'))
    data_differences = eval_shots_against_sql(
        shotlist=shotlist,
        mdsplus_data=mdsplus_data,
        sql_data=sql_data,
        data_columns=[data_column],
        expected_failure_columns=expected_failure_columns,  # we use xfail instead of manually expecting for column failures
        fail_quick=fail_quick,
    )
    if not fail_quick:
        expected_failure = any(
            data_difference.expect_failure for data_difference in data_differences
        )
        if expected_failure:
            pytest.xfail(
                reason="matches expected data failures"
            )  # stops execution of test
        else:
            assert all(
                not data_difference.failed for data_difference in data_differences
            ), get_failure_statistics_string(data_differences, data_column=data_column)


def test_other_values(
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    sql_data: Dict[int, pd.DataFrame],
    data_columns: List[str],
    expected_failure_columns: List[str],
    fail_quick: bool,
):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """

    mdsplus_columns = set().union(*(df.columns for df in mdsplus_data.values()))
    sql_columns = set().union(*(df.columns for df in sql_data.values()))

    test_columns = mdsplus_columns.intersection(sql_columns).difference(data_columns)

    data_differences = eval_shots_against_sql(
        shotlist=shotlist,
        mdsplus_data=mdsplus_data,
        sql_data=sql_data,
        data_columns=test_columns,
        fail_quick=fail_quick,
        expected_failure_columns=expected_failure_columns,
    )

    if not fail_quick:
        matches_expected = all(
            data_difference.matches_expected for data_difference in data_differences
        )
        expected_failure = any(
            data_difference.expect_failure for data_difference in data_differences
        )
        if matches_expected and expected_failure:
            pytest.xfail(
                reason="matches expected data failures"
            )  # stops execution of test
        else:
            assert all(
                not data_difference.failed for data_difference in data_differences
            ), get_failure_statistics_string(data_differences)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fail-slow",
        action="store_true",
        help="Get summary of column failures, for specified column(s)",
    )
    parser.add_argument(
        "--data-column",
        type=str.lower,
        default=None,
        help="Data column to test, use all data columns if not specified",
    )

    parser.add_argument(
        "--shot-id",
        type=int,
        action="store",
        default=None,
        help="Shot number to test, uses the default shot list if not specified",
    )

    args = parser.parse_args()

    fail_quick = not args.fail_slow
    data_columns = [args.data_column] if args.data_column else None
    tokamak = resolve_tokamak_from_environment()

    if args.shot_id is None:
        shotlist = get_tokamak_test_shotlist(tokamak)
    else:
        shotlist = [args.shot_id]

    expected_failure_columns = get_tokamak_test_expected_failure_columns(tokamak)

    data_differences = eval_against_sql(
        tokamak=tokamak,
        shotlist=shotlist,
        expected_failure_columns=expected_failure_columns,
        fail_quick=fail_quick,
        test_columns=data_columns,
    )

    print(get_failure_statistics_string(data_differences))
