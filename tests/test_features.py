#!/usr/bin/env python3

import logging
import os

import pandas as pd
import pytest

from disruption_py.handlers.handler import Handler
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.utils.environment_vars import temporary_env_vars
from disruption_py.utils.exceptions import TokamakNotSupportedError

TEST_SETTINGS = {
    "default_fast": {},
    "sqlcache_full": {"existing_data_request": "sql"},
    "timebase_full": {
        "existing_data_request": "sql",
        "use_existing_data_timebase": True,
        "efit_tree_name": "analysis",
    },
    "columns_full": {
        "existing_data_request": "sql",
        "run_columns": ["v_loop", "q95"],
        "run_tags": [],
        "only_requested_columns": True,
    },
    "flattop_full": {
        "efit_tree_name": "analysis",
        "signal_domain": "flattop",
        "run_tags": [],
        "run_methods": ["_get_ip_parameters"],
    },
    "rampup_fast": {
        "cmod": {
            "efit_tree_name": "analysis",
            "signal_domain": "rampup_and_flattop",
        }
    },
    "logging_full": {
        "log_settings": {
            "log_to_file": True,
            "file_log_level": logging.WARNING,
            "log_file_write_mode": "a",
            "log_to_console": False,
            "console_log_level": logging.WARNING,
            "use_custom_logging": False,
        },
        "efit_tree_name": "efit18",
        "signal_domain": "flattop",
    },
}


@pytest.mark.parametrize("shot_settings_key", TEST_SETTINGS.keys())
def test_features_serial(
    handler: Handler, tokamak, shotlist, shot_settings_key, test_file_path_f
):
    if "GITHUB_ACTIONS" in os.environ and "_fast" not in shot_settings_key:
        pytest.skip("fast execution")

    with temporary_env_vars([("FILE_LOG_PATH", test_file_path_f(".log"))]):
        try:
            test_setting = ShotSettings.from_dict(
                TEST_SETTINGS[shot_settings_key], tokamak=tokamak
            )
        except TokamakNotSupportedError:
            pytest.skip(f"not tested for tokamak {tokamak.value}")

        results = handler.get_shots_data(
            shot_ids_request=shotlist,
            shot_settings=test_setting,
            output_type_request=[
                "list",
                "dataframe",
                test_file_path_f(".csv"),
                test_file_path_f(".hdf5"),
            ],
            num_processes=1,
        )
    list_output, df_output, csv_processed, hdf_processed = results
    assert isinstance(list_output, list)
    assert isinstance(df_output, pd.DataFrame)
    assert csv_processed == hdf_processed == len(shotlist)
    assert os.path.exists(test_file_path_f(".csv"))
    assert os.path.exists(test_file_path_f(".hdf5"))
    assert (
        os.path.exists(test_file_path_f(".log"))
        or not test_setting.log_settings.log_to_file
    )


def test_features_parallel(handler: Handler, tokamak, shotlist, test_file_path_f):
    with temporary_env_vars([("FILE_LOG_PATH", test_file_path_f(".log"))]):
        try:
            test_setting = ShotSettings.from_dict(
                TEST_SETTINGS["default_fast"], tokamak=tokamak
            )
        except TokamakNotSupportedError:
            pytest.skip(f"not tested for tokamak {tokamak.value}")

        results = handler.get_shots_data(
            shot_ids_request=shotlist,
            shot_settings=test_setting,
            output_type_request=[
                "list",
                "dataframe",
                test_file_path_f(".csv"),
                test_file_path_f(".hdf5"),
            ],
            num_processes=2,
        )
    list_output, df_output, csv_processed, hdf_processed = results
    assert isinstance(list_output, list)
    assert isinstance(df_output, pd.DataFrame)
    assert csv_processed == hdf_processed == len(shotlist)
    assert os.path.exists(test_file_path_f(".csv"))
    assert os.path.exists(test_file_path_f(".hdf5"))
    assert (
        os.path.exists(test_file_path_f(".log"))
        or not test_setting.log_settings.log_to_file
    )
