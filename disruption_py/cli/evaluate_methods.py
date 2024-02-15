import argparse
from contextlib import contextmanager
from typing import Dict, List
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import logging
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.shot_ids_request import ShotIdsRequestParams, shot_ids_request_runner
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.utils.constants import TIME_CONST
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment
from disruption_py.utils.math_utils import matlab_gradient_1d_vectorized

CMOD_TEST_SHOTS = [
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

TIME_EPSILON = 0.05 # Tolerance for taking the difference between two times [s]
IP_EPSILON = 1e5    # Tolerance for taking the difference between two ip values [A]

VAL_TOLERANCE = 0.01   # Tolerance for comparing values between MDSplus and SQL
MATCH_FRACTION = 0.95   # Fraction of signals that must match between MDSplus and SQL

def get_mdsplus_data(handler, shot_list):
    shot_settings = ShotSettings(
        efit_tree_name="efit18",
        set_times_request="efit",
        log_settings=LogSettings(
            console_log_level=logging.ERROR
        )
    )
    return handler.get_shots_data(
        shot_ids_request=shot_list, 
        shot_settings=shot_settings,
        output_type_request="dict",
    )
    
def get_sql_data(handler, mdsplus_data : Dict, shot_list):
    shot_data = {}
    for shot_id in shot_list:
        times = mdsplus_data[shot_id]['time']
        sql_data =  handler.database.get_shots_data([shot_id])
        shot_data[shot_id] = pd.merge_asof(times.to_frame(), sql_data, on='time', direction='nearest', tolerance=TIME_CONST)
    return shot_data


def test_data_match(sql_shot_df : pd.DataFrame, mdsplus_shot_df : pd.DataFrame, data_column : str):
    
    sql_column_data = sql_shot_df[data_column].astype(np.float64)
    mds_column_data = mdsplus_shot_df[data_column].astype(np.float64)
    
    # copare data for numeric differences
    relative_difference = np.where(
        sql_column_data != 0, 
        np.abs((mds_column_data - sql_column_data) / sql_column_data), 
        np.where(mds_column_data != 0, np.inf, np.nan)
    )
    numeric_anomalies_mask = np.greater(relative_difference, VAL_TOLERANCE)
    
    # compare data for nan differences
    sql_is_nan_ = pd.isnull(sql_column_data)
    mdsplus_is_nan = pd.isnull(mds_column_data)
    nan_anomalies_mask = (sql_is_nan_ != mdsplus_is_nan)
    
    anomalies = np.argwhere(numeric_anomalies_mask | nan_anomalies_mask)
    
    return not (len(anomalies) / len(relative_difference) > 1 - MATCH_FRACTION)

def evaluate_cmod_accuracy(shot_list : List = None):
    """
    Evaluate the accuracy of CMod methods.
    
    Prints a short report on the methods that have suceeded and failed.
    Success criteria is having more that 95% of results within 1% of known results. 
    """
    print("Evaluating accuracy of CMod methods...")
    if shot_list is None or len(shot_list) == 0:
        shot_list = CMOD_TEST_SHOTS
    else:
        shot_list = [int(shot_id) for shot_id in shot_list]
    
    @contextmanager
    def monkey_patch_numpy_gradient():
        original_function = np.gradient
        np.gradient = matlab_gradient_1d_vectorized
        try:
            yield
        finally:
            np.gradient = original_function
    
    with monkey_patch_numpy_gradient():
        return _evaluate_cmod_accuracy(shot_list)
    
def _evaluate_cmod_accuracy(shot_list : List = None):
    cmod_handler = CModHandler()
    print("Getting data from MDSplus")
    mdsplus_data = get_mdsplus_data(cmod_handler, shot_list)
    print("Getting data from sql table")
    sql_data = get_sql_data(cmod_handler, mdsplus_data, shot_list)
    
    success_columns = set()
    unknown_columns = set()
    failure_columns = set()

    for shot_id in shot_list:
        
        mdsplus_shot_df : pd.DataFrame = mdsplus_data[shot_id]
        sql_shot_df : pd.DataFrame = sql_data[shot_id]
        
        for data_column in mdsplus_shot_df.columns:
            if data_column not in sql_shot_df.columns or not is_numeric_dtype(mdsplus_shot_df[data_column]):
                unknown_columns.add(data_column)
                continue
                
            if test_data_match(sql_shot_df, mdsplus_shot_df, data_column):
                success_columns.add(data_column)
            else:
                failure_columns.add(data_column)
    
    success_columns = success_columns.difference(failure_columns)
    unknown_columns = unknown_columns.difference(success_columns).difference(failure_columns)   
    print(f"Successful Columns (failure criteria not met for any shot): {success_columns}") 
    print(f"Columns with a failure: {failure_columns}") 
    print(f"Columns that lacked testing data: {unknown_columns}") 
    return success_columns, unknown_columns, failure_columns

def main(args):
    """
    Test the accuracy of DisruptionPy on a Tokamak
    
    Prints a short report on the methods that have suceeded and failed.
    Success criteria is having more that 95% of results within 1% of known results. 
    """
    print(f"Welcome to the test accuracy script for disruption_py")
    
    tokamak_string: str = input('Which tokamak would you like to evaluate methods for? ("cmod" for cmod, "d3d" for d3d, leave blank to auto-detect): ') 
    while True:
        if tokamak_string == "":
            tokamak = get_tokamak_from_environment()
        else:
            tokamak = map_string_to_enum(tokamak_string, Tokamak, should_raise=False)
        
        if tokamak is None:
            tokamak_string = input('Invalid input, please enter "cmod" for cmod, "d3d" for d3d, or leave blank to auto-detect: ')
        else:
            break   
    
    # Get shot ids
    if args.shotlist is None:
        print()
        print('What shot ids would you like to test? (Please enter shot ids as comma seperated list and/or on seperate lines)')
        print('To finish please leave a blank line (if no shots are entered a short default shot list will be used)') 
        all_shot_ids = []
        while True:
            shot_ids: str = input()
            if shot_ids == "":
                break
            
            shot_ids_to_add = [feature.strip() for feature in shot_ids.split(",")]
            for shot_id in shot_ids_to_add:
                try:
                    all_shot_ids.append(int(shot_id))
                except ValueError:
                    print(f"Shot id {shot_id} does not appear to be an integer, skipping")
    else:
        shot_ids_request_params = ShotIdsRequestParams(
            database=None,
            tokamak=tokamak,
            logger=logging.getLogger("test_accuracy_logger")
        )
        all_shot_ids = shot_ids_request_runner(args.shotlist, shot_ids_request_params)
    
    if tokamak == Tokamak.CMOD:
        evaluate_cmod_accuracy(all_shot_ids)
    else:
        print("Sorry, this tokamak is not currently supported.")
        
        
def get_parser():
    parser = argparse.ArgumentParser(description='Evaluate the accuracy of DisruptionPy methods on a Tokamak.')
    parser.add_argument('--shotlist', type=str, help='Path to file specifying a shotlist, leave blank for interactive mode', default=None)
    return parser