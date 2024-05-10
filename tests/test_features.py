import logging
import pytest
import os
import pandas as pd

from disruption_py.handlers.handler import Handler
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.shot_settings import ShotSettings

FAST_SHOT_COUNT = 1

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
    "rampup_fast": ShotSettings(
        efit_tree_name="analysis",
        signal_domain="rampup_and_flattop",
    ),
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


@pytest.fixture(scope="module")
def testing_shotlist(shotlist):
    if "GITHUB_ACTIONS" in os.environ:
        return shotlist[:FAST_SHOT_COUNT]
    return shotlist


@pytest.fixture(scope="module")
def shot_settings_keys():
    if "GITHUB_ACTIONS" in os.environ:
        return [k for k in TEST_SETTINGS if k.endswith("_fast")]
    return TEST_SETTINGS.keys()


@pytest.mark.skipif(
    os.path.exists("/fusion/projects/disruption_warning"), reason="on DIII-D"
)
@pytest.mark.parametrize("shot_settings_key", TEST_SETTINGS.keys())
def test_features_serial(
    handler : Handler, testing_shotlist, shot_settings_key, shot_settings_keys
):
    if shot_settings_key not in shot_settings_keys:
        pytest.skip("fast execution")
    results = handler.get_shots_data(
        shot_ids_request=testing_shotlist,
        shot_settings=TEST_SETTINGS[shot_settings_key],
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
    assert csv_processed == hdf_processed == len(testing_shotlist)


@pytest.mark.skipif(
    os.path.exists("/fusion/projects/disruption_warning"), reason="on DIII-D"
)
def test_features_parallel(handler : Handler, testing_shotlist):
    results = handler.get_shots_data(
        shot_ids_request=testing_shotlist,
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
    assert csv_processed == hdf_processed == len(testing_shotlist)


@pytest.fixture(scope="session", autouse=True)
def cleanup_after_tests(request):
    yield
    for ext in ["log", "csv", "hdf5"]:
        delete_file_path = f"{__file__}.{ext}"
        if os.path.exists(delete_file_path):
            os.remove(delete_file_path)
