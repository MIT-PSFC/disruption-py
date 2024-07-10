#!/usr/bin/env python3

import os

import pandas as pd
import pytest

from disruption_py.workflow import get_shots_data
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.machine.tokamak import is_tokamak_indexed

TEST_SETTINGS = {
    "default_fast": {},
    "sqlcache_full": {"cache_setting": "sql"},
    "timebase_full": {
        "cache_setting": "sql",
        "use_cache_setting_timebase": True,
        "efit_tree_name": "analysis",
    },
    "columns_full": {
        "cache_setting": "sql",
        "run_columns": ["v_loop", "q95"],
        "run_tags": [],
        "only_requested_columns": True,
    },
    "flattop_full": {
        "efit_tree_name": "analysis",
        "domain_setting": "flattop",
        "run_tags": [],
        "run_methods": ["_get_ip_parameters"],
    },
    "rampup_fast": {
        "cmod": {
            "efit_tree_name": "analysis",
            "domain_setting": "rampup_and_flattop",
        }
    },
}


@pytest.mark.parametrize("retrieval_settings_key", TEST_SETTINGS.keys())
def test_features_serial(tokamak, shotlist, retrieval_settings_key, test_file_path_f):
    if "GITHUB_ACTIONS" in os.environ and "_fast" not in retrieval_settings_key:
        pytest.skip("fast execution")

    test_setting = TEST_SETTINGS[retrieval_settings_key]
    if is_tokamak_indexed(test_setting):
        if tokamak.value not in test_setting:
            pytest.skip(f"not tested for tokamak {tokamak.value}")

    test_setting = RetrievalSettings.from_dict(test_setting, tokamak=tokamak)

    results = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=test_setting,
        output_setting=[
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


def test_features_parallel(tokamak, shotlist, test_file_path_f):
    test_setting = RetrievalSettings.from_dict(
        TEST_SETTINGS["default_fast"], tokamak=tokamak
    )

    results = get_shots_data(
        shotlist_setting=shotlist,
        retrieval_settings=test_setting,
        output_setting=[
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
