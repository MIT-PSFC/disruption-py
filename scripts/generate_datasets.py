import json
import argparse
import sys
import logging
from disruption_py.handlers.cmod_handler import CModHandler

from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.settings import LogSettings, ShotSettings

try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

import pandas as pd 

from disruption_py.utils.math_utils import generate_id
import disruption_py.data
from disruption_py.ml.preprocessing import *


def main(args):
    if args.log:
        log_settigs = LogSettings(log_file_path=fr'./output/{args.unique_id}.log', file_log_level=args.log_level*10)
    else:
        log_settigs = LogSettings(console_log_level=args.log_level*10)
    logger = log_settigs.logger()

    feature_cols, derived_feature_cols = parse_feature_cols(args.feature_cols)
    logger.info(f"requested feature columns: {feature_cols}")
    logger.info(f"requested derived feature columns: {derived_feature_cols}")
    
    shot_settings = ShotSettings(
        log_settings=log_settigs, 
        efit_tree_name=args.efit_tree,
        set_times_request=args.timebase_signal,
        run_methods=args.populate_methods,
        run_tags=args.populate_tags,
        run_columns=feature_cols,
        only_requested_columns=args.only_requested_columns,
        output_type_request="list"
    )
        
    if args.shotlist is None:
        shot_ids_request = "paper"
    else:
        shot_ids_request = args.shotlist
    
    tokemak = map_string_to_enum(args.tokamak, Tokamak)
    if tokemak == Tokamak.CMOD:
        handler = CModHandler()
    else:
        raise ValueError("Tokamak Not Supported")
    
    dataset_df = handler.get_shots_data(
        shot_ids_request=shot_ids_request, 
        shot_settings=shot_settings,
        num_processes=args.num_processes
    )
    
    if args.filter:
        dataset_df = filter_dataset_df(dataset_df, exclude_non_disruptive=False,
                                        exclude_black_window=BLACK_WINDOW_THRESHOLD, impute=True)
    
    dataset_df = add_derived_features(dataset_df, derived_feature_cols)
    dataset_df['label'] = create_label(dataset_df['time_until_disrupt'])
    dataset_df.to_csv(args.output_dir +
                      f"whole_df_{args.unique_id}.csv", sep=',', index=False)
    
    if args.produce_train_test:
        X_train, X_test, y_train, y_test = create_dataset(
            dataset_df, ratio=DEFAULT_RATIO)
        logger.info(f"{X_train.shape}, {X_test.shape}, {y_train.shape}, {y_test.shape}")
        logger.info(f"Shots (Train): {len(pd.unique(X_train['shot']))}")
        logger.info(f"Shots (Test): {len(pd.unique(X_test['shot']))}")
        df_train = pd.concat([X_train, y_train], axis=1)
        df_train.to_csv(args.output_dir +
                        f"train_{args.unique_id}.csv", sep=',', index=False)
        
        df_test = pd.concat([X_test, y_test], axis=1)
        df_test.to_csv(args.output_dir +
                    f"test_{args.unique_id}.csv", sep=',', index=False)
        
    with open(args.output_dir + f"generate_datasets_{args.unique_id}.json", "w") as f:
        json.dump(vars(args), f)
    logger.info(f"Unique ID for this run: {args.unique_id}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate DPRF compatible datasets for training and inference. Currently only supports DIII-D data")
    parser.add_argument('--shotlist', type=str,
                        help='Path to file specifying shotlist', default=None)
    parser.add_argument('--tokamak', type=str,
                        help='Tokamak to use for data source. Currently only supports DIII-D and Alcator C-Mod.', default='d3d')
    parser.add_argument('--num_processes', type=int,
                        help='The numberof processes to use for data retrieval.', default=1)
    parser.add_argument(
        '--only_requested_columns', type=bool, help="Whether to only create a datset with the requested columns", default=False)
    parser.add_argument('--feature_cols', type=str,
                        help='Either a file or comma-separated list of desired feature columns', default=None)
    parser.add_argument('--output_dir', type=str,
                        help='Path to generated data.', default=r'./output/')
    parser.add_argument('--timebase_signal', type=str,
                        help='Signal whose timebase will be used as the unifying timebase of the dataset.', default=None)
    parser.add_argument('--efit_tree', type=str,
                        help="Name of efit tree to use for each shot. If left as None, the script will use the get_efit_tree method in database.py.", default=None)
    parser.add_argument('--data_source', type=int, choices=[
                        0, 1, 2, 3], help=r"0: Default to SQL database then MDSPlus.\n1: Default to MDSPlus then SQL database.\n2: SQL database only.\n3: MDSPlus only.", default=2)
    parser.add_argument('--unique_id', type=str,
                        help='Unique identifier for the dataset. Used to name the output files.', default=generate_id())
    parser.add_argument(
        '--log', type=bool, help='By default, generate_datasets will log to commandline but if this argument is true it will log to a file in the output directory', default=False)
    parser.add_argument('--log_level', type=int, choices=[
                        0, 1, 2, 3, 4, 5], help='Notset:0,Debug:1,Info:2,Warning:3,Error:4,Critical:5', default=2)
    parser.add_argument('--label', type=str, choices=[
                        'binary', 'none'], help="Timestep disruption label. Currently only supports binary labels", default='binary')
    # Two argparse arguments. populate_methods which is a list of strings but defaults to None and populate_tags which is a list of strings that also defaults to None
    parser.add_argument('--populate_methods', nargs='*', type=str, default=None)
    parser.add_argument('--populate_tags', nargs='*', type=str, default=None) 
    parser.add_argument(
        '--filter', type=bool, help="Run filter_dataset method on produced dataset. Necessary for generating DPRF datasets", default=True)
    parser.add_argument(
        '--produce_train_test', type=bool, help="Whether to proucee train and test datasets", default=False)
    args = parser.parse_args()
    main(args)
