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
from disruption_py.utils.eval.data_difference import DataDifference
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment
from disruption_py.utils.eval.test_helpers import get_mdsplus_data, get_sql_data_for_mdsplus, test_against_sql
from disruption_py.utils.eval.environment_constants import get_test_handler, get_test_expected_failure_columns, get_test_shot_ids

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
    if data_column in expected_failure_columns:
        pytest.xfail(f"Known failure for {data_column}")
    
    data_differences = DataDifference.test_shots(
        shot_ids=shotlist,
        mdsplus_data=mdsplus_data, 
        sql_data=sql_data, 
        data_columns=[data_column],
        expected_failure_columns=[], # we use xfail instead of manually expecting for column failures
        fail_quick=fail_quick
    )
    
    if not fail_quick:
        assert all((not data_difference.failed) for data_difference in data_differences), DataDifference.get_failure_statistics_string(
        data_differences, data_column=data_column)
    
    
def test_other_values(shotlist : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], tested_data_columns : List[str], expected_failure_columns : List[str], fail_quick : bool):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    
    mdsplus_columns = set().union(*(df.columns for df in mdsplus_data.values()))
    sql_columns = set().union(*(df.columns for df in sql_data.values()))
    
    test_columns = mdsplus_columns.intersection(sql_columns).difference(tested_data_columns)
    
    data_differences = DataDifference.test_shots(
        shot_ids=shotlist, 
        mdsplus_data=mdsplus_data, 
        sql_data=sql_data, 
        data_columns=test_columns,
        fail_quick=fail_quick,
        expected_failure_columns=expected_failure_columns
    )
    
    if not fail_quick:
        if any(data_column in expected_failure_columns for data_column in test_columns):
            pytest.xfail(f"Known failure for at least one column")
        
        assert all((not data_difference.failed) for data_difference in data_differences), DataDifference.get_failure_statistics_string(data_differences)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fail-quick', action='store_true')
    parser.add_argument('--data-column', default=None, help='Data column to test, use all data columns if not specified')
    args = parser.parse_args()
    
    fail_quick = args.fail_quick
    data_columns = [args.data_column] if args.data_column else None
    tokamak = get_tokamak_from_environment()
    
    handler = get_test_handler(tokamak)
    shot_ids = get_test_shot_ids(tokamak)
    expected_failure_columns = get_test_expected_failure_columns(tokamak)
    
    
    test_against_sql(handler=handler, shot_ids=shot_ids, expected_failure_columns=expected_failure_columns, fail_quick=fail_quick, test_columns=data_columns)

