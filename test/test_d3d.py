"""Unit tests for workflows involving get_dataset_df() for obtaining D3D data.

Expects to be run on the MFE workstations.
Expects MDSplus to be installed and configured.
Expects SQL credentials to be configured.
"""
from typing import Dict
from disruption_py.handlers.d3d_handler import D3DHandler
import pytest

import numpy as np
import pandas as pd
import logging
import os
from disruption_py.handlers.d3d_handler import D3DHandler
from disruption_py.settings import ShotSettings, LogSettings
from disruption_py.settings.set_times_request import ListSetTimesRequest 
from disruption_py.utils.constants import TIME_CONST 

USER = os.getenv('USER')

# Shot list used for testing
# Mix of disruptive and non-disruptive shots present in SQL and MDSplus
TEST_SHOTS = [
    161228, # disruptive
    # 161237, # disruptive
    # 166177, # non disruptive 
    # 166253
]

TEST_COLUMNS = [
    'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt',
    'beta_p', 'beta_n', 'li', 'n_equal_1_mode_IRLM', 'z_error', 'v_z',
    'kappa', 'H98', 'q0', 'qstar', 'q95', 'dn_dt', 'radiated_fraction',
    'power_supply_railed', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt',
    'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'dWmhd_dt',
    'dprad_dt', 'p_nbi', 'p_ech', 'p_ohm', 'intentional_disruption',
    'Greenwald_fraction', 'Te_HWHM', 'other_hardware_failure', 'Te_HWHM_RT',
    'v_loop_RT', 'n_e_RT', 'Greenwald_fraction_RT', 'ip_error_RT', 'ip_RT',
    'dipprog_dt_RT', 'Wmhd_RT', 'Wmhd', 'n_equal_1_mode',
    'n_equal_1_normalized', 'Te_width_normalized', 'Te_width_normalized_RT',
    'q95_RT', 'li_RT', 'beta_p_RT', 'oeamp1em', 'oeamp1om', 'oefrq1em',
    'oefrq1om', 'oeamp1e', 'oeamp1o', 'oefrq1e', 'oefrq1o', 'delta',
    'squareness', 'zcur_normalized', 'aminor', 'n1rms_normalized',
    'kappa_area', 'Te_peaking_CVA_RT', 'ne_peaking_CVA_RT',
    'Prad_peaking_CVA_RT', 'Prad_peaking_XDIV_RT', 'H_alpha',
]

KNOWN_FAILURE_COLUMNS = []

# TEST_COLUMNS = list(set(TEST_COLUMNS).difference(KNOWN_FAILURE_COLUMNS))

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

VAL_TOLERANCE = 0.01   # Tolerance for comparing values between MDSplus and SQL
MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL

@pytest.fixture(scope='module')
def d3d_handler():
    return D3DHandler()

@pytest.fixture(scope='module')
def shotlist():
    return TEST_SHOTS

@pytest.fixture(scope='module')
def mdsplus_data(d3d_handler : D3DHandler, shotlist) -> Dict:
    shot_settings = ShotSettings(
        set_times_request="disruption",
        log_settings=LogSettings(
            log_to_console=False,
            log_file_path=f"/cscratch/{USER}/last_test_log.log",
            log_file_write_mode="w",
            file_log_level=logging.DEBUG
        )
    )
    shot_data = d3d_handler.get_shots_data(
        shot_ids_request=shotlist, 
        shot_settings=shot_settings,
        output_type_request="dict",
    )
    return shot_data

@pytest.fixture(scope='module')
def sql_data(d3d_handler : D3DHandler, shotlist, mdsplus_data : Dict):
    shot_data = {}
    for shot_id in shotlist:
        times = mdsplus_data[shot_id]['time']
        sql_data =d3d_handler.database.get_shots_data([shot_id])
        assert len(times) == len(sql_data), f"Shot {shot_id} has {len(times)} rows but SQL has {len(sql_data)} rows"
        
        shot_data[shot_id] = pd.merge_asof(times.to_frame(), sql_data, on='time', direction='nearest', tolerance=TIME_CONST)
        assert len(times) == len(shot_data[shot_id]), f"After projecting timebase shot {shot_id} has {len(times)} rows but SQL has {len(shot_data[shot_id])} rows"
        
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
    

# 41 33, 57 17
# FAILED test/test_d3d.py::test_data_columns[0-False-n_equal_1_mode_IRLM] - AssertionError: Column n_equal_1_mode_IRLM missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-v_z] - AssertionError: Column v_z missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-H98] - AssertionError: Shot 161228 column H98 failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-radiated_fraction] - AssertionError: Column radiated_fraction missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-zcur] - AssertionError: Shot 161228 column zcur failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-v_loop] - AssertionError: Column v_loop missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-p_rad] - AssertionError: Column p_rad missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-dprad_dt] - AssertionError: Column dprad_dt missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-p_nbi] - AssertionError: Column p_nbi missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-p_ech] - AssertionError: Column p_ech missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-p_ohm] - AssertionError: Column p_ohm missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-intentional_disruption] - AssertionError: Column intentional_disruption missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-Te_HWHM] - AssertionError: Column Te_HWHM missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-other_hardware_failure] - AssertionError: Column other_hardware_failure missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-Te_HWHM_RT] - AssertionError: Column Te_HWHM_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-v_loop_RT] - AssertionError: Column v_loop_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-Greenwald_fraction_RT] - AssertionError: Shot 161228 column Greenwald_fraction_RT failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-dipprog_dt_RT] - AssertionError: Shot 161228 column dipprog_dt_RT failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-n_equal_1_mode] - AssertionError: Shot 161228 column n_equal_1_mode failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-n_equal_1_normalized] - AssertionError: Shot 161228 column n_equal_1_normalized failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-Te_width_normalized] - AssertionError: Column Te_width_normalized missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-Te_width_normalized_RT] - AssertionError: Column Te_width_normalized_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-q95_RT] - AssertionError: Shot 161228 column q95_RT failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-li_RT] - AssertionError: Shot 161228 column li_RT failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-beta_p_RT] - AssertionError: Shot 161228 column beta_p_RT failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-oeamp1em] - AssertionError: Column oeamp1em missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oeamp1om] - AssertionError: Column oeamp1om missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oefrq1em] - AssertionError: Column oefrq1em missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oefrq1om] - AssertionError: Column oefrq1om missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oeamp1e] - AssertionError: Column oeamp1e missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oeamp1o] - AssertionError: Column oeamp1o missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oefrq1e] - AssertionError: Column oefrq1e missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-oefrq1o] - AssertionError: Column oefrq1o missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-zcur_normalized] - AssertionError: Shot 161228 column zcur_normalized failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-n1rms_normalized] - AssertionError: Shot 161228 column n1rms_normalized failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-kappa_area] - AssertionError: Shot 161228 column kappa_area failed for arrays:
# FAILED test/test_d3d.py::test_data_columns[0-False-Te_peaking_CVA_RT] - AssertionError: Column Te_peaking_CVA_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-ne_peaking_CVA_RT] - AssertionError: Column ne_peaking_CVA_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-Prad_peaking_CVA_RT] - AssertionError: Column Prad_peaking_CVA_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-Prad_peaking_XDIV_RT] - AssertionError: Column Prad_peaking_XDIV_RT missing from MDSPlus for shot 161228
# FAILED test/test_d3d.py::test_data_columns[0-False-H_alpha] - AssertionError: Shot 161228 column H_alpha failed for arrays: