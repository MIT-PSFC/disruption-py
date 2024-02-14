
import logging
import pytest
import os
import pandas as pd

from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.shot_settings import ShotSettings

TEST_SHOTS = [
    1150805012,   # Flattop Disruption
    1150805013,     # No Disruption
    1150805014,     # No Disruption
    1150805015,     # Rampdown Disruption
    1150805016,     # Rampdown Disruption
]

@pytest.fixture(scope='module')
def cmod_handler():
    return CModHandler()

@pytest.fixture(scope='module')
def shot_settings_list():
    return [
        ShotSettings(),
          ShotSettings(
            existing_data_request="sql",
        ),
        ShotSettings(
            existing_data_request="sql",
               use_existing_data_timebase = True,
              efit_tree_name="analysis",
        ),
        ShotSettings(
            existing_data_request="sql",
               run_columns=["v_loop", "q95"],
            run_tags=[],
            only_requested_columns=True,
        ),
         ShotSettings(
            efit_tree_name="analysis",
            signal_domain = "flattop",
               run_tags=[],
               run_methods=["_get_ip_parameters"]
        ),
          ShotSettings(
            efit_tree_name="analysis",
            signal_domain = "rampup_and_flattop",
        ),
        ShotSettings(
            # logging
            log_settings=LogSettings(
                log_file_path="test/temp/test_log.log",
                file_log_level=logging.WARNING,
                log_file_write_mode="a",
                log_to_console=False,
                console_log_level=logging.WARNING,
                use_custom_logging=False,
            ),
            
            # data settings
            efit_tree_name="efit18",
            signal_domain = "flattop",
        )
  
    ]


@pytest.mark.parametrize('multiprocessing', [True, False])
@pytest.mark.parametrize('shot_settings_index', list(range(7)))
def test_features(cmod_handler, shot_settings_list, shot_settings_index, multiprocessing):
    list_output, df_output, num_processed, num_processed = cmod_handler.get_shots_data(
        shot_ids_request=TEST_SHOTS,
        shot_settings=shot_settings_list[shot_settings_index],
        output_type_request=["list", "dataframe", "test/temp/test_output.csv", "test/temp/test_output.hdf5"],
        num_processes=4 if multiprocessing else 1
    )
    assert isinstance(list_output, list)
    assert isinstance(df_output, pd.DataFrame)
    
@pytest.fixture(scope="session", autouse=True)
def cleanup_after_tests(request):
    yield
    # Teardown code: delete files after all tests are done
    delete_file_paths = ["test/temp/test_log.log", "test/temp/test_output.csv", "test/temp/test_output.hdf5"]
    for delete_file_path in delete_file_paths:
        if os.path.exists(delete_file_path):
            os.remove(delete_file_path)