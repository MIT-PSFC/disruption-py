import json
import argparse
import sys
import logging
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

import pandas as pd 

from disruption_py.utils import generate_id
import disruption_py.data
from disruption_py.ml.preprocessing import *


def main(args):
    if args.log:
        ch = logging.FileHandler(fr'./output/{args.unique_id}.log')
    else:
        ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(args.log_level*10)
    # log_format = '%(asctime)s | %(levelname)s: %(message)s'
    # ch.setFormatter(logging.Formatter(log_format))
    LOGGER.addHandler(ch)
    LOGGER.setLevel(args.log_level*10)
    print(LOGGER)
    if args.shotlist is None:
        with importlib_resources.path(disruption_py.data, "paper_shotlist.txt") as p:
            args.shotlist = str(p)
    shot_ids = pd.read_csv(
        args.shotlist, header=None).iloc[:, 0].values.tolist()
    print(shot_ids)
    feature_cols, derived_feature_cols = parse_feature_cols(args.feature_cols)
    print(feature_cols)
    dataset_df = get_dataset_df(args.data_source, cols=feature_cols +
                                REQUIRED_COLS, efit_tree=args.efit_tree, shot_ids=shot_ids, tokamak=args.tokamak, timebase_signal=args.timebase_signal, label=args.label, populate=args.populate)
    dataset_df = add_derived_features(dataset_df, derived_feature_cols)
    if args.filter is True:
        print(args.filter)
        dataset_df = filter_dataset_df(dataset_df, exclude_non_disruptive=False,
                                       exclude_black_window=BLACK_WINDOW_THRESHOLD, impute=True)
    X_train, X_test, y_train, y_test = create_dataset(
        dataset_df, ratio=DEFAULT_RATIO)
    dataset_df.to_csv(args.output_dir +
                      f"whole_df_{args.unique_id}.csv", sep=',', index=False)
    # df_train_val = pd.concat([X_train, y_train], axis=1)
    # X_train, X_val, y_train, y_val = create_dataset(df_train_val, ratio=.25)
    print(X_train.shape, X_test.shape,
          y_train.shape, y_test.shape)
    df_train = pd.concat([X_train, y_train], axis=1)
    df_train.to_csv(args.output_dir +
                    f"train_{args.unique_id}.csv", sep=',', index=False)
    # df_val = pd.concat([X_val, y_val], axis=1)
    # df_val.to_csv(args.output_dir +
    # f"val_{args.unique_id}.csv", sep=",", index=False)
    df_test = pd.concat([X_test, y_test], axis=1)
    df_test.to_csv(args.output_dir +
                   f"test_{args.unique_id}.csv", sep=',', index=False)
    with open(args.output_dir + f"generate_datasets_{args.unique_id}.json", "w") as f:
        json.dump(vars(args), f)
    print(f"Unique ID for this run: {args.unique_id}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate DPRF compatible datasets for training and inference. Currently only supports DIII-D data")
    parser.add_argument('--shotlist', type=str,
                        help='Path to file specifying shotlist', default=None)
    parser.add_argument('--tokamak', type=str,
                        help='Tokamak to use for data source. Currently only supports DIII-D and Alcator C-Mod.', default='d3d')
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
    parser.add_argument('--populate', type=str, choices=[
                        'default', 'l_mode'], help="Options for population at instantiation.", default='default')
    parser.add_argument(
        '--filter', type=int, help="Run filter_dataset method on produced dataset. Necessary for generating DPRF datasets", default=1)
    args = parser.parse_args()
    main(args)
