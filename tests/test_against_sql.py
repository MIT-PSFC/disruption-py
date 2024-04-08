"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
from typing import Dict, List
import pytest

import numpy as np
import pandas as pd
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.utils.eval.data_difference import DataDifference
from disruption_py.utils.eval.test_helpers import get_mdsplus_data
from disruption_py.utils.eval.test_helpers import get_sql_data_for_mdsplus
from disruption_py.utils.mappings.tokamak_helpers import get_handler_from_environment 

# Shot list used for testing
# Mix of disruptive and non-disruptive shots present in SQL and MDSplus
TEST_SHOTS = [
    1150805012,   # Flattop Disruption
    1150805013,     # No Disruption
    1150805014,     # No Disruption
    1150805015,     # Rampdown Disruption
    1150805016,     # Rampdown Disruption
    1150805017,     # Rampdown Disruption
    1150805019,     # Rampdown Disruption
    1150805020,     # Rampdown Disruption
    1150805021,     # Rampdown Disruption
    1150805022      # Flattop Disruption
]

TEST_COLUMNS = [
    'I_efc', 'sxr', 'time_until_disrupt', 'beta_n', 'beta_p', 'kappa', 'li',
    'upper_gap', 'lower_gap', 'q0', 'qstar', 'q95', 'v_loop_efit', 'Wmhd',
    'ssep', 'n_over_ncrit', 'tritop', 'tribot', 'a_minor', 'rmagx', 'chisq',
    'dbetap_dt', 'dli_dt', 'dWmhd_dt', 'V_surf', 'kappa_area', 'Te_width',
    'ne_peaking', 'Te_peaking', 'pressure_peaking', 'n_e', 'dn_dt',
    'Greenwald_fraction', 'n_equal_1_mode', 'n_equal_1_normalized',
    'n_equal_1_phase', 'BT', 'prad_peaking', 'v_0', 'ip', 'dip_dt',
    'dip_smoothed', 'ip_prog', 'dipprog_dt', 'ip_error', 'z_error',
    'z_prog', 'zcur', 'v_z', 'z_times_v_z', 'p_oh', 'v_loop', 'p_rad',
    'dprad_dt', 'p_lh', 'p_icrf', 'p_input', 'radiated_fraction', 'time',
    'shot', 'commit_hash'
]

KNOWN_FAILURE_COLUMNS = [
    'lower_gap', 'upper_gap', 'ssep', 'dipprog_dt', 'n_over_ncrit', # constant factor scaling error
    'ip_error' # constant error
]

# TEST_COLUMNS = list(set(TEST_COLUMNS).difference(KNOWN_FAILURE_COLUMNS))

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

@pytest.fixture(scope='module')
def handler_fixture():
    return get_handler_from_environment()

@pytest.fixture(scope='module')
def shotlist_fixture():
    return TEST_SHOTS

@pytest.fixture(scope='module')
def mdsplus_data_fixture(handler_fixture : CModHandler, shotlist_fixture : List[int]) -> Dict[int, pd.DataFrame]:
    return get_mdsplus_data(handler_fixture, shotlist_fixture)

@pytest.fixture(scope='module')
def sql_data_fixture(handler_fixture : CModHandler, shotlist_fixture : List[int], mdsplus_data_fixture : Dict[int, pd.DataFrame]) -> Dict[int, pd.DataFrame]:
    return get_sql_data_for_mdsplus(handler_fixture, shotlist_fixture, mdsplus_data_fixture)

@pytest.mark.parametrize("data_column", TEST_COLUMNS)
def test_data_columns(shotlist_fixture : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], data_column, verbose_output):
    data_differences = DataDifference.test_shots(shotlist_fixture, mdsplus_data, sql_data, [data_column])
    
    assert all((not data_difference.failed) for data_difference in data_differences), DataDifference.get_failure_statistics_string(
        data_differences, verbose_output, data_column=data_column)
    
def test_other_values(shotlist_fixture : List[int], mdsplus_data : Dict[int, pd.DataFrame], sql_data : Dict[int, pd.DataFrame], verbose_output):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    
    mdsplus_columns = set().union(*(df.columns for df in mdsplus_data.values()))
    sql_columns = set().union(*(df.columns for df in sql_data.values()))
    
    test_columns = mdsplus_columns.intersection(sql_columns).difference(TEST_COLUMNS)
        
    data_differences = DataDifference.test_shots(shotlist_fixture, mdsplus_data, sql_data, test_columns)
    
    assert all((not data_difference.failed) for data_difference in data_differences), DataDifference.get_failure_statistics_string(data_differences)


   
    

