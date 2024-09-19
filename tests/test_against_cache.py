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
    eval_against_cache,
    eval_shots_against_cache,
    get_cached_from_fresh,
    get_fresh_data,
)
from tests.utils.factory import (
    get_tokamak_test_expected_failure_columns,
    get_tokamak_test_shotlist,
)
from tests.utils.pytest_helper import extract_param, save_to_csv


@pytest.fixture(scope="module", name="fresh_data")
def fresh_data_fixture(
    tokamak: Tokamak,
    shotlist: List[int],
    module_file_path_f,
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    fresh_data = get_fresh_data(
        tokamak=tokamak,
        shotlist=shotlist,
        log_file_path=module_file_path_f(".log"),
        test_columns=extract_param(pytestconfig),
    )
    save_to_csv(
        data=fresh_data, module_file_path_f=module_file_path_f, data_source_name="fresh"
    )
    return fresh_data


@pytest.fixture(scope="module", name="cache_data")
def cache_data_fixture(
    tokamak: Tokamak,
    shotlist: List[int],
    fresh_data: Dict[int, pd.DataFrame],
    module_file_path_f,
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    cache_data = get_cached_from_fresh(
        tokamak=tokamak,
        shotlist=shotlist,
        fresh_data=fresh_data,
        test_columns=extract_param(pytestconfig),
    )
    save_to_csv(
        data=cache_data, module_file_path_f=module_file_path_f, data_source_name="cache"
    )
    return cache_data


def test_data_columns(
    shotlist: List[int],
    fresh_data: Dict[int, pd.DataFrame],
    cache_data: Dict[int, pd.DataFrame],
    data_column,
    expected_failure_columns: List[str],
):
    """
    Test that the data columns are the same between fresh and cached sources across
    specified data columns.

    Data column is parameterized in pytest_generate_tests.
    """
    eval_shots_against_cache(
        shotlist=shotlist,
        fresh_data=fresh_data,
        cache_data=cache_data,
        data_columns=[data_column],
        expected_failure_columns=expected_failure_columns,  # we use xfail instead of manually expecting for column failures
    )


def main():
    """
    main function called by command-line invocation.
    """

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

    data_differences = eval_against_cache(
        tokamak=tokamak,
        shotlist=shotlist,
        expected_failure_columns=expected_failure_columns,
        test_columns=data_columns,
    )

    columns = {dd.data_column for dd in data_differences}
    print(
        f"Python tests complete. Checked {len(shotlist)} shots with {len(columns)} columns."
    )


if __name__ == "__main__":
    main()
