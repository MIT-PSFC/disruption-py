"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
import argparse
from typing import Dict, List
import pytest

import pandas as pd
from disruption_py.handlers.cmod_handler import Handler
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment, get_test_handler, get_test_expected_failure_columns, get_test_shot_ids
from disruption_py.utils.eval.eval_against_sql import eval_shots_against_sql, get_failure_statistics_string, get_mdsplus_data, get_sql_data_for_mdsplus, eval_against_sql

@pytest.fixture(scope='module')
def mdsplus_data(handler : Handler, shotlist : List[int]) -> Dict[int, pd.DataFrame]:
    return get_mdsplus_data(handler, shotlist)

@pytest.fixture(scope='module')
def sql_data(handler : Handler, shotlist : List[int], mdsplus_data : Dict[int, pd.DataFrame]) -> Dict[int, pd.DataFrame]:
    return get_sql_data_for_mdsplus(handler, shotlist, mdsplus_data)

def test_data_columns(shotlist : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], data_column, expected_failure_columns : List[str], fail_quick : bool):
    """
    Test that the data columns are the same between MDSplus and SQL across specified data columns.
    
    Data column is parameterized in pytest_generate_tests.
    """
    # if data_column in expected_failure_columns:
    #     request.node.add_marker(pytest.mark.xfail(reason='column expected failure'))
    data_differences = eval_shots_against_sql(
        shot_ids=shotlist,
        mdsplus_data=mdsplus_data, 
        sql_data=sql_data, 
        data_columns=[data_column],
        expected_failure_columns=expected_failure_columns, # we use xfail instead of manually expecting for column failures
        fail_quick=fail_quick
    )
    if not fail_quick:
        matches_expected = all(data_difference.matches_expected for data_difference in data_differences)
        expected_failure = any(data_difference.expect_failure for data_difference in data_differences)
        if matches_expected and expected_failure:
            pytest.xfail(reason='matches expected data failures') # stops execution of test
        else:
            assert all(not data_difference.failed for data_difference in data_differences), get_failure_statistics_string(
                data_differences, data_column=data_column)
    
    
def test_other_values(shotlist : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], data_columns : List[str], expected_failure_columns : List[str], fail_quick : bool):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    
    mdsplus_columns = set().union(*(df.columns for df in mdsplus_data.values()))
    sql_columns = set().union(*(df.columns for df in sql_data.values()))
    
    test_columns = mdsplus_columns.intersection(sql_columns).difference(data_columns)
    
    data_differences = eval_shots_against_sql(
        shot_ids=shotlist, 
        mdsplus_data=mdsplus_data, 
        sql_data=sql_data, 
        data_columns=test_columns,
        fail_quick=fail_quick,
        expected_failure_columns=expected_failure_columns
    )
    
    if not fail_quick:
        matches_expected = all(data_difference.matches_expected for data_difference in data_differences)
        expected_failure = any(data_difference.expect_failure for data_difference in data_differences)
        if matches_expected and expected_failure:
            pytest.xfail(reason='matches expected data failures') # stops execution of test
        else:
            assert all(not data_difference.failed for data_difference in data_differences), get_failure_statistics_string(data_differences)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fail-slow', action='store_true', help="Get summary of column failures, for specified column(s)")
    parser.add_argument('--data-column', default=None, help='Data column to test, use all data columns if not specified')
    args = parser.parse_args()
    
    fail_quick = not args.fail_slow
    data_columns = [args.data_column] if args.data_column else None
    tokamak = get_tokamak_from_environment()
    
    handler = get_test_handler(tokamak)
    shot_ids = get_test_shot_ids(tokamak)
    expected_failure_columns = get_test_expected_failure_columns(tokamak)
    
    data_differences = eval_against_sql(
        handler=handler, 
        shot_ids=shot_ids, 
        expected_failure_columns=expected_failure_columns, 
        fail_quick=fail_quick, 
        test_columns=data_columns
    )
    
    print(get_failure_statistics_string(data_differences))

