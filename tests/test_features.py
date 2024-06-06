#!/usr/bin/env python3

import logging
import os

import pandas as pd
import pytest

from disruption_py.handlers.handler import Handler
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.shot_settings import ShotSettings

TEST_SETTINGS = {
    "default_fast": ShotSettings(),
    "sqlcache_full": ShotSettings(existing_data_request="sql"),
    "timebase_full": ShotSettings(
        existing_data_request="sql",
        use_existing_data_timebase=True,
        efit_tree_name="analysis",
    ),
    "columns_full": ShotSettings(
        existing_data_request="sql",
        run_columns=["v_loop", "q95"],
        run_tags=[],
        only_requested_columns=True,
    ),
    "flattop_full": ShotSettings(
        efit_tree_name="analysis",
        signal_domain="flattop",
        run_tags=[],
        run_methods=["_get_ip_parameters"],
    ),
    "rampup_fast": {
        "cmod": ShotSettings(
            efit_tree_name="analysis",
            signal_domain="rampup_and_flattop",
        )
    },
    "logging_full": ShotSettings(
        log_settings=LogSettings(
            log_file_path=f"{__file__}.log",
            file_log_level=logging.WARNING,
            log_file_write_mode="a",
            log_to_console=False,
            console_log_level=logging.WARNING,
            use_custom_logging=False,
        ),
        efit_tree_name="efit18",
        signal_domain="flattop",
    ),
}


@pytest.mark.parametrize("shot_settings_key", TEST_SETTINGS.keys())
def test_features_serial(handler: Handler, tokamak, shotlist, shot_settings_key):
    if "GITHUB_ACTIONS" in os.environ and "_fast" not in shot_settings_key:
        pytest.skip("fast execution")

    test_setting = TEST_SETTINGS[shot_settings_key]
    if isinstance(test_setting, dict):
        if tokamak.value in test_setting:
            test_setting = test_setting[tokamak.value]
        else:
            pytest.skip(f"not tested for tokamak {tokamak.value}")

    results = handler.get_shots_data(
        shot_ids_request=shotlist,
        shot_settings=test_setting,
        output_type_request=[
            "list",
            "dataframe",
            f"{__file__}.csv",
            f"{__file__}.hdf5",
        ],
        num_processes=1,
    )
    list_output, df_output, csv_processed, hdf_processed = results
    assert isinstance(list_output, list)
    assert isinstance(df_output, pd.DataFrame)
    assert csv_processed == hdf_processed == len(shotlist)


def test_features_parallel(handler: Handler, shotlist):
    results = handler.get_shots_data(
        shot_ids_request=shotlist,
        shot_settings=TEST_SETTINGS["default_fast"],
        output_type_request=[
            "list",
            "dataframe",
            f"{__file__}.csv",
            f"{__file__}.hdf5",
        ],
        num_processes=2,
    )
    list_output, df_output, csv_processed, hdf_processed = results
    assert isinstance(list_output, list)
    assert isinstance(df_output, pd.DataFrame)
    assert csv_processed == hdf_processed == len(shotlist)


@pytest.fixture(scope="session", autouse=True)
def cleanup_after_tests(request):
    yield
    for ext in ["log", "csv", "hdf5"]:
        delete_file_path = f"{__file__}.{ext}"
        if os.path.exists(delete_file_path):
            os.remove(delete_file_path)
