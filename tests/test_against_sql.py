#!/usr/bin/env python3

"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""

import argparse
from typing import Dict, List

import pandas as pd
import pytest

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from tests.utils.eval_against_sql import (
    eval_against_sql,
    eval_shots_against_sql,
    get_mdsplus_data,
    get_sql_data_for_mdsplus,
)
from tests.utils.factory import (
    get_tokamak_test_expected_failure_columns,
    get_tokamak_test_shotlist,
)

from tests.utils.pytest_helper import extract_param, save_to_csv


@pytest.fixture(scope="module")
def mdsplus_data(
    tokamak: Tokamak,
    shotlist: List[int],
    module_file_path_f,
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    mds = get_mdsplus_data(
        tokamak=tokamak,
        shotlist=shotlist,
        log_file_path=module_file_path_f(".log"),
        test_columns=extract_param(pytestconfig),
    )
    save_to_csv(data=mds, module_file_path_f=module_file_path_f, data_source_name="mds")
    return mds


@pytest.fixture(scope="module")
def sql_data(
    tokamak: Tokamak,
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    module_file_path_f,
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    sql = get_sql_data_for_mdsplus(
        tokamak=tokamak,
        shotlist=shotlist,
        mdsplus_data=mdsplus_data,
        test_columns=extract_param(pytestconfig),
    )
    save_to_csv(data=sql, module_file_path_f=module_file_path_f, data_source_name="sql")
    return sql


def test_data_columns(
    shotlist: List[int],
    mdsplus_data: Dict[int, pd.DataFrame],
    sql_data: Dict[int, pd.DataFrame],
    data_column,
    expected_failure_columns: List[str],
):
    """
    Test that the data columns are the same between MDSplus and SQL across specified data columns.

    Data column is parameterized in pytest_generate_tests.
    """
    eval_shots_against_sql(
        shotlist=shotlist,
        mdsplus_data=mdsplus_data,
        sql_data=sql_data,
        data_columns=[data_column],
        expected_failure_columns=expected_failure_columns,  # we use xfail instead of manually expecting for column failures
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        test_columns=data_columns,
    )

    columns = {dd.data_column for dd in data_differences}
    print(
        f"Python tests complete. Checked {len(shotlist)} shots with {len(columns)} columns."
    )
