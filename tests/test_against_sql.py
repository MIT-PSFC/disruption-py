"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
from typing import Dict, List
import pytest

import numpy as np
import pandas as pd
from disruption_py.handlers.cmod_handler import Handler
from disruption_py.utils.constants import CMOD_TEST_SHOTS
from disruption_py.utils.eval.data_difference import DataDifference
from disruption_py.utils.eval.test_helpers import get_mdsplus_data
from disruption_py.utils.eval.test_helpers import get_sql_data_for_mdsplus
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_handler_from_environment 

KNOWN_FAILURE_COLUMNS = [
    'lower_gap', 'upper_gap', 'ssep', 'dipprog_dt', 'n_over_ncrit', # constant factor scaling error
    'ip_error' # constant error
]

# TEST_COLUMNS = list(set(TEST_COLUMNS).difference(KNOWN_FAILURE_COLUMNS))

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]
    
@pytest.fixture(scope='module')
def mdsplus_data_fixture(handler : Handler, shotlist : List[int]) -> Dict[int, pd.DataFrame]:
    return get_mdsplus_data(handler, shotlist)

@pytest.fixture(scope='module')
def sql_data_fixture(handler : Handler, shotlist : List[int], mdsplus_data_fixture : Dict[int, pd.DataFrame]) -> Dict[int, pd.DataFrame]:
    return get_sql_data_for_mdsplus(handler, shotlist, mdsplus_data_fixture)

def test_data_columns(shotlist : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], data_column, verbose_output):
    """
    Test that the data columns are the same between MDSplus and SQL across specified data columns.
    
    Data column is parameterized in pytest_generate_tests.
    """
    data_differences = DataDifference.test_shots(shotlist, mdsplus_data, sql_data, [data_column])
    
    assert all((not data_difference.failed) for data_difference in data_differences), DataDifference.get_failure_statistics_string(
        data_differences, verbose_output, data_column=data_column)
    
def test_other_values(shotlist : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], data_columns : List[str]):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    
    mdsplus_columns = set().union(*(df.columns for df in mdsplus_data.values()))
    sql_columns = set().union(*(df.columns for df in sql_data.values()))
    
    test_columns = mdsplus_columns.intersection(sql_columns).difference(data_columns)
        
    data_differences = DataDifference.test_shots(shotlist, mdsplus_data, sql_data, test_columns)
    
    assert all((not data_difference.failed) for data_difference in data_differences), DataDifference.get_failure_statistics_string(data_differences)


   
    

