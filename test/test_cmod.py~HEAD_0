"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
import pytest

import numpy as np
import pandas as pd

from disruption_py.shots import CModShot 
from disruption_py.database import create_cmod_handler 

# Shot list used for testing
# Mix of disruptive and non-disruptive shots present in SQL and MDSplus
TEST_SHOTS = [1150805012,   # Flattop Disruption
            1150805013,     # No Disruption
            1150805014,     # No Disruption
            1150805015,     # Rampdown Disruption
            1150805016,     # Rampdown Disruption
            1150805017,     # Rampdown Disruption
            1150805019,     # Rampdown Disruption
            1150805020,     # Rampdown Disruption
            1150805021,     # Rampdown Disruption
            1150805022]     # Flattop Disruption
TEST_COLS = ['beta_n', 'beta_p', 'kappa', 'li', 'upper_gap', 'lower_gap', 'q0',
       'qstar', 'q95', 'v_loop_efit', 'Wmhd', 'ssep', 'n_over_ncrit', 'v_surf',
       'R0', 'tritop', 'tribot', 'a_minor', 'chisq', 'dbetap_dt', 'dli_dt',
       'dWmhd_dt', 'H98', 'Te_width', 'n_e', 'dn_dt', 'Greenwald_fraction',
       'I_efc', 'ip', 'dip_dt', 'dip_smoothed', 'ip_prog', 'dipprog_dt',
       'ip_error', 'kappa_area', 'n_equal_1_mode', 'n_equal_1_normalized',
      'n_equal_1_phase', 'BT', 'p_oh', 'v_loop', 'p_rad', 'dprad_dt', 'p_lh',
       'p_icrf', 'p_input', 'radiated_fraction', 'v_0', 'sxr',
       'time_until_disrupt']
TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL

@pytest.fixture(scope='module')
def cmod():
    return create_cmod_handler()

def check_sql_columns(cmod):
    """
    Ensure there is no error when getting all columns from the SQL table for several shots.
    """
    sql_df = cmod.get_shot_data(TEST_SHOTS, cols='*')
    assert sql_df is not None

@pytest.fixture(scope='module')
def shotlists(cmod):
    expected_shots = []
    test_shots = []
    for shot_id in TEST_SHOTS:
        test_shot = CModShot(shot_id, disruption_time = cmod.get_disruption_time(shot_id))
        expected_shot = cmod.get_shot(shot_id) 
        test_shot.data.sort_values('time', inplace=True)
        expected_shot.data.sort_values('time', inplace=True)
        merged_into_expected = pd.merge_asof(expected_shot.data,test_shot.data, on='time', direction='nearest', tolerance=0.001)
        merged_into_shot = pd.merge_asof(test_shot.data,expected_shot.data, on='time', direction='nearest', tolerance=0.001)
        merged_into_expected = merged_into_expected.dropna(axis=0, subset=['shot_x','shot_y'])
        merged_into_shot = merged_into_shot.dropna(axis=0, subset=['shot_x','shot_y'])
        matched_times_expected = merged_into_shot['time']
        matched_times_shot = merged_into_expected['time']
        test_shot.data = test_shot.data[test_shot.data['time'].isin(matched_times_shot)].reset_index(drop=True)
        expected_shot.data = expected_shot.data[expected_shot.data['time'].isin(matched_times_expected)].reset_index(drop=True)
        assert len(test_shot.data) == len(expected_shot.data), f"Shot {test_shot} has {len(test_shot.data)} rows but SQL has {len(expected_shot.data)} rows"
        expected_shots.append(expected_shot)
        test_shots.append(test_shot)
    return test_shots, expected_shots

def test_get_all_columns(shotlists):
    test_shots, _ = shotlists
    for col in TEST_COLS:
        for shot in test_shots:
            assert col in shot.data.columns, f"Shot {shot} is missing {col} from MDSplus"

@pytest.mark.parametrize("col", TEST_COLS)
def test_parameter_calculations(shotlists, col):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    test_shots, expected_shots = shotlists
    for shot,expected_shot in zip(test_shots, expected_shots):
        if col not in expected_shot.data.columns:
            continue
        # check if the col of the shot is all nan
        if np.isnan(shot.data[col]).all():
            if np.isnan(expected_shot.data[col]).all():
                continue
        # Compare percentage diff between MDSplus and SQL. In the case that the SQL value is 0, inf should be the diff if the MDSplus value is non-zero and nan if the MDSplus value is 0
        diff = np.where(expected_shot.data[col] != 0, 
                        np.abs((shot.data[col] - expected_shot.data[col]) / expected_shot.data[col]), 
                        np.where(shot.data[col] != 0, np.inf, np.nan))                       
        anomalies = np.argwhere(diff > 1.e-2)
        assert len(anomalies) / len(diff) < MATCH_FRACTION, f"Shot {shot} has {len(anomalies)} anomalies for {col} out of {len(diff)}"

# Other tests for MDSplus
# TODO: Refactor these tests

@pytest.mark.parametrize("shot", TEST_SHOTS[:1])
def test_flattop_times(cmod, shot):
    """Ensure that the flattop time matches the change in programmed ip in the SQL dataset.
    """
    sql_df = cmod.get_shot(shot).data[['time', 'dipprog_dt']]
    # Find the first time where dipprog_dt is zero
    sql_flattop_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[0]
    # Find the last time where dipprog_dt is zero
    sql_flattop_end_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[-1]

    flattop_df = CModShot(shot, disruption_time = cmod.get_disruption_time(shot), timebase_signal='flattop').data[['time', 'dipprog_dt']]
    
    # Find the first time in the flattop signal
    mds_flattop_time = flattop_df['time'].iloc[0]
    # Find the last time in the flattop signal
    mds_flattop_end_time = flattop_df['time'].iloc[-1]

    assert mds_flattop_time == pytest.approx(sql_flattop_time, abs=TIME_EPSILON)
    assert mds_flattop_end_time == pytest.approx(sql_flattop_end_time, abs=TIME_EPSILON)