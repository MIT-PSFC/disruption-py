"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
import pytest

import numpy as np
import pandas as pd
import logging
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings import ShotSettings, LogSettings
from disruption_py.settings.set_times_request import ListSetTimesRequest 
from disruption_py.utils.constants import TIME_CONST 

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

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

VAL_TOLERANCE = 0.01   # Tolerance for comparing values between MDSplus and SQL
MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL

@pytest.fixture(scope='module')
def cmod():
    return CModHandler()

def get_mdsplus_data(cmod_handler: CModHandler, shot_id):
    shot_settings = ShotSettings(
        efit_tree_name="efit18",
        set_times_request="efit",
        log_settings=LogSettings(
            log_to_console=False,
            log_file_path="test/last_log.log",
            file_log_level=logging.DEBUG
        )
    )
    shot_data = cmod_handler.get_shots_data(shot_id, shot_settings=shot_settings)
    return shot_data[0]

def get_sql_data(cmod_handler: CModHandler, shot_id, times):
    sql_data =cmod_handler.database.get_shot_data([shot_id])
    sql_data["sql_time"] = sql_data["time"]
    return pd.merge_asof(times.to_frame(), sql_data, on='time', direction='nearest', tolerance=TIME_CONST)

@pytest.fixture(scope='module')
def shotlists(cmod):
    expected_shots = []
    test_shots = []
    for shot_id in TEST_SHOTS:
        test_shot_data = get_mdsplus_data(cmod, shot_id)
        expected_shot_data = get_sql_data(cmod, shot_id, test_shot_data['time'])
        assert len(test_shot_data) == len(expected_shot_data), f"Shot {shot_id} has {len(test_shot_data)} rows but SQL has {len(expected_shot_data)} rows"
        expected_shots.append(expected_shot_data)
        test_shots.append(test_shot_data)
    return test_shots, expected_shots

SKIPPABLE_COLUMNS = [
    'lower_gap', 'upper_gap', 'ssep', 'n_over_ncrit', 'dipprog_dt', # constant factor scaling error
    'ip_error' # unknown error
]

# lower_gap, upper_gap, ssep is times 100 in sql table
            
@pytest.mark.parametrize("fail_early", [True, False])
def test_all_sql_values(shotlists, fail_early):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    test_shots, expected_shots = shotlists
    successful_shot_cols = []
    failed_shot_cols = []
    for shot_id, test_shot_data, expected_shot_data in zip(TEST_SHOTS, test_shots, expected_shots):
        for col in expected_shot_data.columns:
            # uncomment to skip columns that are recognized to be broken
            # if col in SKIPPABLE_COLUMNS:
            #     continue
            
            if col not in test_shot_data.columns:
                print(f"Shot {shot_id} is missing {col} from MDSplus source")
                continue
            
            # check if the col of the shot is all nan
            if test_shot_data[col].isna().all() and expected_shot_data[col].isna().all():
                continue
            
            # Compare percentage diff between MDSplus and SQL. In the case that the SQL value is 0, inf should be the diff if the MDSplus value is non-zero and nan if the MDSplus value is 0
            diff = np.where(expected_shot_data[col] != 0, 
                            np.abs((test_shot_data[col] - expected_shot_data[col]) / expected_shot_data[col]), 
                            np.where(test_shot_data[col] != 0, np.inf, np.nan))                       
            anomalies = np.argwhere(diff > VAL_TOLERANCE)
            if len(anomalies) / len(diff) > 1 - MATCH_FRACTION:
                failed_shot_cols.append((shot_id, col, len(anomalies), len(diff), f"{(len(anomalies)/len(diff)*100):.2f}%"))
                if fail_early:
                    indexes = np.arange(len(diff)) # anomalies.flatten()
                    anomaly_differences = diff[indexes]
                    test_shot_data_differences = test_shot_data[col].iloc[indexes]
                    expected_shot_data_differences = expected_shot_data[col].iloc[indexes]
                    anomaly = np.where(diff > VAL_TOLERANCE, 1, 0)[indexes]
                    difference_df = pd.DataFrame({'Test': test_shot_data_differences, 'expected': expected_shot_data_differences, 'difference': anomaly_differences, 'anomaly': anomaly})
                    # difference_df.to_csv(f"test/cmod_failed_values_{shot_id}_{col}.csv")
                    pd.options.display.max_rows = None
                    raise AssertionError(
                        f"Shot {shot_id} condition failed for {col}. Arrays:\n{difference_df}"
                    )
            else:
                successful_shot_cols.append((shot_id, col, len(anomalies), len(diff), f"{(len(anomalies)/len(diff)*100):.2f}%"))
                    
    if len(failed_shot_cols) > 0:
        success_cols = {successful_shot_col[1] for successful_shot_col in successful_shot_cols}
        failure_cols = {failed_shot_col[1] for failed_shot_col in failed_shot_cols}
        print(f"Succeeded on columns: {success_cols - failure_cols}")
        print(f"Failed on columns: {failure_cols}")

        col_success_percent = len(successful_shot_cols) / (len(successful_shot_cols) + len(failed_shot_cols)) * 100
        col_success_string = f"Succeeded on {col_success_percent:.2f}% of columns\n"
        
        total_anomalies = sum([x[2] for x in failed_shot_cols]) + sum([x[2] for x in successful_shot_cols])
        total_values = sum([x[3] for x in failed_shot_cols]) + sum([x[3] for x in successful_shot_cols])
        entry_success_percent = (total_values - total_anomalies) / total_values * 100
        entry_success_string = f"Succeeded on {entry_success_percent:.2f}% of entries\n"
        
        sorted_failed_with_shot_cols = sorted(failed_shot_cols, key=lambda x: (x[2] / x[3]), reverse=True)
        list_of_anomalies_string = f"Failed with anomalies (shot_id, col, num_anomalies, num_differences):\n{sorted_failed_with_shot_cols}"
        raise AssertionError("\n" + col_success_string + entry_success_string + list_of_anomalies_string)
        
# Other tests for MDSplus
# TODO: Refactor these tests

@pytest.mark.parametrize("shot_id", TEST_SHOTS[:1])
def test_flattop_times(cmod, shot_id):
    """
    Ensure that the flattop time matches the change in programmed ip in the SQL dataset.
    """
    sql_df = get_sql_data(cmod, shot_id)[['time', 'dipprog_dt']]
    # Find the first time where dipprog_dt is zero
    sql_flattop_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[0]
    # Find the last time where dipprog_dt is zero
    sql_flattop_end_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[-1]

    shot_settings = ShotSettings(
        set_times_request='efit',
        signal_domain='flattop'
    )
    
    flattop_df = cmod.get_shots_data(shot_id, shot_settings=shot_settings)[0][['time', 'dipprog_dt']]
    
    # Find the first time in the flattop signal
    mds_flattop_time = flattop_df['time'].iloc[0]
    # Find the last time in the flattop signal
    mds_flattop_end_time = flattop_df['time'].iloc[-1]

    assert mds_flattop_time == pytest.approx(sql_flattop_time, abs=TIME_EPSILON)
    assert mds_flattop_end_time == pytest.approx(sql_flattop_end_time, abs=TIME_EPSILON)


# test specific column error in detail
# def  test_derrivatives(shotlists):
#     test_shots, expected_shots = shotlists
#     for shot_id, test_shot_data, expected_shot_data in zip(TEST_SHOTS, test_shots, expected_shots):
#         difference_dfs = []
#         for col in ['beta_p', 'dbetap_dt']:
#             diff = np.where(expected_shot_data[col] != 0, 
#                             np.abs((test_shot_data[col] - expected_shot_data[col]) / expected_shot_data[col]), 
#                             np.where(test_shot_data[col] != 0, np.inf, np.nan))                       
#             indexes = np.arange(len(diff)) # anomalies.flatten()
#             anomaly_differences = diff[indexes]
#             test_shot_data_differences = test_shot_data[col].iloc[indexes]
#             expected_shot_data_differences = expected_shot_data[col].iloc[indexes]
#             anomaly = np.where(diff > VAL_TOLERANCE, 1, 0)[indexes]
#             difference_df = pd.DataFrame({f'Test_{col}': test_shot_data_differences, f'Expected_{col}': expected_shot_data_differences, f'Difference_{col}': anomaly_differences, f'Anomaly_{col}': anomaly})
#             difference_dfs.append(difference_df)
            
#         total_difference_df = pd.concat(difference_dfs, axis=1)
#         # total_difference_df['evals'] = np.gradient(total_difference_df['Expected_beta_p'], np.round(expected_shot_data['time'], 3), edge_order=1)
#         # total_difference_df['evals_diff'] = np.where(total_difference_df['Expected_beta_p'] != 0, 
#         #                     np.abs((total_difference_df['evals'] - total_difference_df['Expected_dbetap_dt']) / total_difference_df['Expected_dbetap_dt']), 
#         #                     np.where(total_difference_df['evals'] != 0, np.inf, np.nan))
#         # total_difference_df['evals_anomaly'] = np.where(total_difference_df['evals_diff'] > VAL_TOLERANCE, 1, 0)
#         total_difference_df['beta_p_diff'] = total_difference_df['Expected_beta_p'].diff()
#         total_difference_df['time'] = expected_shot_data['time']
#         total_difference_df['time_diff'] = expected_shot_data['time'].diff()

#         total_difference_df.to_csv(f"test/cmod_failed_values_{shot_id}_dbetap_dt.csv")
#         raise AssertionError(
#             f"Shot {shot_id} condition failed. Arrays:\n{total_difference_df}"
#         )


# test specified columns

# TEST_COLS = ['beta_n', 'beta_p', 'kappa', 'li', 'upper_gap', 'lower_gap', 'q0',
#        'qstar', 'q95', 'v_loop_efit', 'Wmhd', 'ssep', 'n_over_ncrit', 'V_surf',
#        'R0', 'tritop', 'tribot', 'a_minor', 'chisq', 'dbetap_dt', 'dli_dt',
#        'dWmhd_dt', 'H98', 'Te_width', 'n_e', 'dn_dt', 'Greenwald_fraction',
#        'I_efc', 'ip', 'dip_dt', 'dip_smoothed', 'ip_prog', 'dipprog_dt',
#        'ip_error', 'kappa_area', 'n_equal_1_mode', 'n_equal_1_normalized',
#       'n_equal_1_phase', 'BT', 'p_oh', 'v_loop', 'p_rad', 'dprad_dt', 'p_lh',
#        'p_icrf', 'p_input', 'radiated_fraction', 'v_0', 'sxr',
#        'time_until_disrupt']

# def test_get_all_columns(shotlists):
#     test_shots, _ = shotlists
#     for col in TEST_COLS:
#         for shot_data in test_shots:
#             assert col in shot_data.columns, f"Shot {shot_data} is missing {col} from MDSplus"

# @pytest.mark.parametrize("col", TEST_COLS)
# def test_parameter_calculations(shotlists, col):
#     """
#     Ensure that all parameters are calculated correctly in the MDSplus shot object.
#     """
#     test_shots, expected_shots = shotlists
#     for shot_id, test_shot_data, expected_shot_data in zip(TEST_SHOTS, test_shots, expected_shots):
#         if col not in expected_shot_data.columns:
#             continue
#         # check if the col of the shot is all nan
#         if test_shot_data[col].isna().all() and expected_shot_data[col].isna().all():
#             continue
        
#         # Compare percentage diff between MDSplus and SQL. In the case that the SQL value is 0, inf should be the diff if the MDSplus value is non-zero and nan if the MDSplus value is 0
#         diff = np.where(expected_shot_data != 0, 
#                         np.abs((test_shot_data[col] - expected_shot_data[col]) / expected_shot_data[col]), 
#                         np.where(test_shot_data[col] != 0, np.inf, np.nan))                       
#         anomalies = np.argwhere(diff > 1.e-2)
#         assert len(anomalies) / len(diff) < MATCH_FRACTION, f"Shot {shot_id} has {len(anomalies)} anomalies for {col} out of {len(diff)}"