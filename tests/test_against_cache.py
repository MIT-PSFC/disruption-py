#!/usr/bin/env python3

"""
Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

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
from tests.utils.pytest_helper import extract_param


@pytest.fixture(scope="module", name="fresh_data")
def fresh_data_fixture(
    tokamak: Tokamak,
    shotlist: List[int],
    test_folder_m: str,
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    """
    Fixture to retrieve fresh data for the specified tokamak and shotlist.

    Parameters
    ----------
    tokamak : Tokamak
        The tokamak object used to retrieve data.
    shotlist : List[int]
        The list of shot identifiers to retrieve data for.
    test_folder_m : str
        Output folder.
    pytestconfig : Config
        The pytest configuration object.

    Returns
    -------
    Dict[int, pd.DataFrame]
        A dictionary mapping shot identifiers to their corresponding fresh DataFrames.
    """
    return get_fresh_data(
        tokamak=tokamak,
        shotlist=shotlist,
        folder=test_folder_m,
        test_columns=extract_param(pytestconfig),
    )


@pytest.fixture(scope="module", name="cache_data")
def cache_data_fixture(
    tokamak: Tokamak,
    shotlist: List[int],
    fresh_data: Dict[int, pd.DataFrame],
    pytestconfig,
) -> Dict[int, pd.DataFrame]:
    """
    Fixture to retrieve cached data based on fresh data for the specified tokamak and shotlist.

    Parameters
    ----------
    tokamak : Tokamak
        The tokamak object used to retrieve data.
    shotlist : List[int]
        The list of shot identifiers to retrieve data for.
    fresh_data : Dict[int, pd.DataFrame]
        The fresh data retrieved for the specified shotlist.
    pytestconfig : Config
        The pytest configuration object.

    Returns
    -------
    Dict[int, pd.DataFrame]
        A dictionary mapping shot identifiers to their corresponding cached DataFrames.
    """
    return get_cached_from_fresh(
        tokamak=tokamak,
        shotlist=shotlist,
        fresh_data=fresh_data,
        test_columns=extract_param(pytestconfig),
    )


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
        expected_failure_columns=expected_failure_columns,
    )


def main():
    """
    main function called by command-line invocation.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--data-column",
        type=str.lower,
        default=None,
        help="Data column to test, use all data columns if not specified",
    )

    parser.add_argument(
        "-s",
        "--shot-id",
        type=int,
        action="store",
        default=None,
        help="Shot number to test, uses the default shot list if not specified",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        type=str.lower,
        action="store",
        default="INFO",
        help="Console log level",
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
        console_log_level=args.log_level,
    )

    columns = {dd.data_column for dd in data_differences}
    print(
        f"Python tests complete. Checked {len(shotlist)} shots with {len(columns)} columns."
    )


if __name__ == "__main__":
    main()
