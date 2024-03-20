"""Unit tests for workflows involving get_dataset_df() for obtaining CMOD data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
from typing import Dict
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

KNOWN_NUMERIC_FAILURE_COLUMNS = [
    'lower_gap', 'upper_gap', 'ssep', 'dipprog_dt', 'n_over_ncrit', # constant factor scaling error
    'ip_error' # constant error
]

# TEST_COLUMNS = list(set(TEST_COLUMNS).difference(KNOWN_NUMERIC_FAILURE_COLUMNS))

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

VAL_TOLERANCE = 0.01   # Tolerance for comparing values between MDSplus and SQL
MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL

@pytest.fixture(scope='module')
def cmod_handler():
    return CModHandler()

@pytest.fixture(scope='module')
def shotlist():
    return TEST_SHOTS

@pytest.fixture(scope='module')
def mdsplus_data(cmod_handler : CModHandler, shotlist) -> Dict:
    shot_settings = ShotSettings(
        efit_tree_name="efit18",
        set_times_request="efit",
        log_settings=LogSettings(
            log_to_console=False,
            log_file_path="test/last_log.log",
            log_file_write_mode="w",
            file_log_level=logging.DEBUG
        )
    )
    shot_data = cmod_handler.get_shots_data(
        shot_ids_request=shotlist, 
        shot_settings=shot_settings,
        output_type_request="dict",
    )
    return shot_data

@pytest.fixture(scope='module')
def sql_data(cmod_handler : CModHandler, shotlist, mdsplus_data : Dict):
    shot_data = {}
    for shot_id in shotlist:
        times = mdsplus_data[shot_id]['time']
        sql_data =cmod_handler.database.get_shots_data([shot_id])
        assert len(times) == len(sql_data), f"Shot {shot_id} has {len(times)} rows but SQL has {len(sql_data)} rows"
        shot_data[shot_id] = pd.merge_asof(times.to_frame(), sql_data, on='time', direction='nearest', tolerance=TIME_CONST)
        
    return shot_data

@pytest.mark.parametrize("data_column", TEST_COLUMNS)
def test_data_columns(shotlist, mdsplus_data : Dict, sql_data : Dict, data_column, verbose_output, fail_slow):
    anomaly_ratios = []
    for shot_id in shotlist:
        mdsplus_shot_data, sql_shot_data = mdsplus_data[shot_id], sql_data[shot_id]
        
        if data_column not in mdsplus_shot_data:
            raise AssertionError(f"Column {data_column} missing from MDSPlus for shot {shot_id}")
        
        if data_column not in sql_shot_data:
            raise AssertionError(f"Column {data_column} missing from SQL for shot {shot_id}")
        
        anomaly_ratio = evaluate_differences(
            shot_id=shot_id,
            sql_shot_data=sql_shot_data,
            mdsplus_shot_data=mdsplus_shot_data,
            data_column=data_column,
            verbose_output=verbose_output,
            fail_slow=fail_slow,
        )
        anomaly_ratios.append(anomaly_ratio)
    
    if any(anomaly_ratio['failed'] for anomaly_ratio in anomaly_ratios):
        raise AssertionError(get_failure_statistics_string(anomaly_ratios, verbose_output, data_column=data_column))
    
def test_other_values(shotlist, mdsplus_data : Dict, sql_data : Dict, verbose_output, fail_slow):
    """
    Ensure that all parameters are calculated correctly in the MDSplus shot object.
    """
    anomaly_ratios = []
    for shot_id in shotlist:
        mdsplus_shot_data, sql_shot_data = mdsplus_data[shot_id], sql_data[shot_id]
        mdsplus_unmatched_cols = list(mdsplus_shot_data.columns.difference(sql_shot_data.columns))
        print(f"Shot {shot_id} is missing {mdsplus_unmatched_cols} from SQL source")
        sql_unmatched_cols = list(sql_shot_data.columns.difference(mdsplus_shot_data.columns))
        print(f"Shot {shot_id} is missing {sql_unmatched_cols} from MDSPlus source")

        for data_column in sql_shot_data.columns.intersection(mdsplus_shot_data.columns):
            
            if data_column in TEST_COLUMNS:
                continue
            
            # check if the col of the shot is all nan
            if mdsplus_shot_data[data_column].isna().all() and sql_shot_data[data_column].isna().all():
                continue
            
            anomaly_ratio = evaluate_differences(
                shot_id=shot_id,
                sql_shot_data=sql_shot_data,
                mdsplus_shot_data=mdsplus_shot_data,
                data_column=data_column,
                verbose_output=verbose_output,
                fail_slow=fail_slow,
            )
            anomaly_ratios.append(anomaly_ratio)

                    
    if any(anomaly_ratio['failed'] for anomaly_ratio in anomaly_ratios):
        
        raise AssertionError(get_failure_statistics_string(anomaly_ratios, verbose_output))

def evaluate_differences(shot_id, sql_shot_data, mdsplus_shot_data, data_column, verbose_output, fail_slow):
    # Compare percentage diff between MDSplus and SQL. In the case that the SQL value is 0, inf should be the diff if the MDSplus value is non-zero and nan if the MDSplus value is 0
    relative_difference = np.where(
        sql_shot_data[data_column] != 0, 
        np.abs((mdsplus_shot_data[data_column] - sql_shot_data[data_column]) / sql_shot_data[data_column]), 
        np.where(mdsplus_shot_data[data_column] != 0, np.inf, np.nan)
    )
    numeric_anomalies_mask = (relative_difference > VAL_TOLERANCE)
    
    sql_is_nan_ = pd.isnull(sql_shot_data[data_column])
    mdsplus_is_nan = pd.isnull(mdsplus_shot_data[data_column])
    nan_anomalies_mask = (sql_is_nan_ != mdsplus_is_nan)
    
    anomalies = np.argwhere(numeric_anomalies_mask | nan_anomalies_mask)
    
    if len(anomalies) / len(relative_difference) > 1 - MATCH_FRACTION:
        if fail_slow:
            failed = True
        else:
            indexes = np.arange(len(relative_difference)) if verbose_output else anomalies.flatten()
            anomaly = np.where(relative_difference > VAL_TOLERANCE, 1, 0)[indexes]
            difference_df = pd.DataFrame({
                'MDSplus Data': mdsplus_shot_data[data_column].iloc[indexes], 
                'Reference Data (SQL)': sql_shot_data[data_column].iloc[indexes], 
                'Relative difference': relative_difference[indexes], 
                'Anomaly': anomaly
            })
            pd.options.display.max_rows = None if verbose_output else 10
            raise AssertionError(f"Shot {shot_id} column {data_column} failed for arrays:\n{difference_df}")
    else:
        failed = False
    
    anomaly_ratio = {
        'failed': failed,
        'shot': shot_id,
        'data_column': data_column,
        'anomalies': len(anomalies), 
        'timebase_length': len(relative_difference),
        'failure_percentage' : f"{len(anomalies) / len(relative_difference*100):.2f}",
    }
    return anomaly_ratio

def get_failure_statistics_string(anomaly_ratios, verbose_output, data_column=None):
    anomaly_ratio_by_column = {}
    for anomaly_ratio in anomaly_ratios:
        anomaly_ratio_by_column.setdefault(anomaly_ratio['data_column'], []).append(anomaly_ratio)
    
    failure_strings = {}
    for ratio_data_column, column_anomaly_ratios in anomaly_ratio_by_column.items():
        failures = [anomaly_ratio['shot'] for anomaly_ratio in column_anomaly_ratios if anomaly_ratio['failed']]
        failed = len(failures) > 0
        if not verbose_output and not failed:
            continue
        successes = [anomaly_ratio['shot'] for anomaly_ratio in column_anomaly_ratios if not anomaly_ratio['failed']]
        anomaly_count = sum([anomaly_ratio['anomalies'] for anomaly_ratio in column_anomaly_ratios])
        timebase_count = sum([anomaly_ratio['timebase_length'] for anomaly_ratio in column_anomaly_ratios])
        failure_strings[ratio_data_column] = f"""
        Column {ratio_data_column} {"FAILED" if failed else "SUCCEEDED"} for shot {anomaly_ratio['shot']}
        Failed for {len(failures)} shots: {failures}
        Succeeded for {len(successes)} shots: {successes}
        Total Entry Failure Rate: {anomaly_count / timebase_count * 100:.2f}%
        """
    
    if data_column is not None:
        return failure_strings.get(data_column, "")
    else:
        return '\n'.join(failure_strings.values())
    

# Other tests for MDSplus
# TODO: Refactor these tests

# @pytest.mark.parametrize("shot_id", TEST_SHOTS[:1])
# def test_flattop_times(cmod, shot_id):
#     """
#     Ensure that the flattop time matches the change in programmed ip in the SQL dataset.
#     """
#     sql_df = get_sql_data(cmod, shot_id)[['time', 'dipprog_dt']]
#     # Find the first time where dipprog_dt is zero
#     sql_flattop_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[0]
#     # Find the last time where dipprog_dt is zero
#     sql_flattop_end_time = sql_df['time'].loc[sql_df['dipprog_dt'] == 0].iloc[-1]

#     shot_settings = ShotSettings(
#         set_times_request='efit',
#         signal_domain='flattop'
#     )
    
#     flattop_df = cmod.get_shots_data(shot_id, shot_settings=shot_settings)[0][['time', 'dipprog_dt']]
    
#     # Find the first time in the flattop signal
#     mds_flattop_time = flattop_df['time'].iloc[0]
#     # Find the last time in the flattop signal
#     mds_flattop_end_time = flattop_df['time'].iloc[-1]

#     assert mds_flattop_time == pytest.approx(sql_flattop_time, abs=TIME_EPSILON)
#     assert mds_flattop_end_time == pytest.approx(sql_flattop_end_time, abs=TIME_EPSILON)


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