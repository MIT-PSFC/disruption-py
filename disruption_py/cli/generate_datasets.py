#!/usr/bin/env python3

import argparse
import json

import pandas as pd

from disruption_py.cli.ml.constants import (
    BLACK_WINDOW_THRESHOLD,
    DEFAULT_COLS,
    DEFAULT_RATIO,
)
from disruption_py.cli.ml.preprocessing import (
    add_derived_features,
    create_dataset,
    create_label,
    filter_dataset_df,
    parse_feature_cols,
)
from disruption_py.core.utils.math import generate_id
from disruption_py.core.utils.misc import without_duplicates
from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data


def main(args):
    if args.log:
        log_settigs = LogSettings(
            log_file_path=rf"./output/{args.unique_id}.log",
            file_log_level=args.log_level * 10,
        )
    else:
        log_settigs = LogSettings(console_log_level=args.log_level * 10)
    logger = log_settigs.logger()

    feature_cols, derived_feature_cols = parse_feature_cols(args.feature_cols)
    feature_cols = without_duplicates(DEFAULT_COLS + feature_cols)
    logger.info(f"requested feature columns: {feature_cols}")
    logger.info(f"requested derived feature columns: {derived_feature_cols}")

    retrieval_settings = RetrievalSettings(
        efit_tree_name=args.efit_tree,
        time_setting=args.timebase_signal,
        run_methods=args.run_methods,
        run_tags=args.run_tags,
        run_columns=feature_cols,
        only_requested_columns=args.only_requested_columns,
        cache_setting="sql" if args.data_source == 0 else None,
    )

    tokamak = resolve_tokamak_from_environment(args)

    if args.shotlist is None:
        if tokamak == Tokamak.D3D:
            shotlist_setting = "d3d_paper_shotlist"
        elif tokamak == Tokamak.CMOD:
            shotlist_setting = "cmod_test"

    dataset_df = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist_setting,
        retrieval_settings=retrieval_settings,
        output_setting="dataframe",
        num_processes=args.num_processes,
        log_settings=log_settigs,
    )

    if args.filter:
        dataset_df = filter_dataset_df(
            dataset_df,
            exclude_non_disruptive=False,
            exclude_black_window=BLACK_WINDOW_THRESHOLD,
            impute=True,
        )

    dataset_df = add_derived_features(dataset_df, derived_feature_cols)
    dataset_df["label"] = create_label(dataset_df["time_until_disrupt"])
    dataset_df.to_csv(
        args.output_dir + f"whole_df_{args.unique_id}.csv", sep=",", index=False
    )

    if args.produce_train_test:
        X_train, X_test, y_train, y_test = create_dataset(
            dataset_df, ratio=DEFAULT_RATIO
        )
        logger.info(f"{X_train.shape}, {X_test.shape}, {y_train.shape}, {y_test.shape}")
        logger.info(f"Shots (Train): {len(pd.unique(X_train['shot']))}")
        logger.info(f"Shots (Test): {len(pd.unique(X_test['shot']))}")
        df_train = pd.concat([X_train, y_train], axis=1)
        df_train.to_csv(
            args.output_dir + f"train_{args.unique_id}.csv", sep=",", index=False
        )

        df_test = pd.concat([X_test, y_test], axis=1)
        df_test.to_csv(
            args.output_dir + f"test_{args.unique_id}.csv", sep=",", index=False
        )
    with open(args.output_dir + f"generate_datasets_{args.unique_id}.json", "w") as f:
        args_dict = vars(args)
        keys_to_remove = ["func", "subcommand"]
        for key in keys_to_remove:
            if key in args_dict:
                del args_dict[key]
        json.dump(args_dict, f)
    logger.info(f"Unique ID for this run: {args.unique_id}")


def get_parser():
    parser = argparse.ArgumentParser(
        description="Generate DPRF compatible datasets for training and inference. Currently only supports CMod data."
    )
    parser.add_argument(
        "--shotlist", type=str, help="Path to file specifying shotlist", default=None
    )
    parser.add_argument(
        "--tokamak",
        type=str,
        help='Tokamak to use for data source. Currently only supports Alcator C-Mod ("cmod" for cmod).',
        default=None,
    )
    parser.add_argument(
        "--num_processes",
        type=int,
        help="The numberof processes to use for data retrieval.",
        default=1,
    )
    parser.add_argument(
        "--only_requested_columns",
        type=bool,
        help="Whether to only create a datset with the requested columns",
        default=False,
    )
    parser.add_argument(
        "--feature_cols",
        type=str,
        help="Either a file or comma-separated list of desired feature columns. Similar to run columns in `RetrievalSettings`",
        default=None,
    )
    parser.add_argument(
        "--output_dir", type=str, help="Path to generated data.", default=r"./output/"
    )
    parser.add_argument(
        "--timebase_signal",
        type=str,
        help="Signal whose timebase will be used as the unifying timebase of the dataset.",
        default="efit",
    )
    parser.add_argument(
        "--efit_tree",
        type=str,
        help="Name of efit tree to use for each shot.",
        default="analysis",
    )
    parser.add_argument(
        "--data_source",
        type=int,
        choices=[0, 1],
        help=r"0: Default to SQL database then MDSPlus.\n1: MDSPlus only.",
        default=1,
    )
    parser.add_argument(
        "--unique_id",
        type=str,
        help="Unique identifier for the dataset. Used to name the output files.",
        default=generate_id(),
    )
    parser.add_argument(
        "--log",
        type=bool,
        help="By default, generate_datasets will log to commandline but if this argument is true it will log to a file in the output directory",
        default=False,
    )
    parser.add_argument(
        "--log_level",
        type=int,
        choices=[0, 1, 2, 3, 4, 5],
        help="Notset:0,Debug:1,Info:2,Warning:3,Error:4,Critical:5",
        default=2,
    )
    parser.add_argument(
        "--label",
        type=str,
        choices=["binary", "none"],
        help="Timestep disruption label. Currently only supports binary labels",
        default="binary",
    )
    # Two argparse arguments. populate_methods which is a list of strings but defaults to None and populate_tags which is a list of strings that also defaults to None
    parser.add_argument("--run_methods", nargs="*", type=str, default=None)
    parser.add_argument("--run_tags", nargs="*", type=str, default=None)
    parser.add_argument(
        "--filter",
        type=bool,
        help="Run filter_dataset method on produced dataset. Necessary for generating DPRF datasets",
        default=True,
    )
    parser.add_argument(
        "--produce_train_test",
        type=bool,
        help="Whether to proucee train and test datasets",
        default=False,
    )
    return parser
