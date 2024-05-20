import logging

from disruption_py.databases.cmod_database import CModDatabase
from disruption_py.handlers import CModHandler
from disruption_py.settings import LogSettings, ShotSettings

handler = CModHandler(
    database_initializer=CModDatabase.default,
    mds_connection_str="alcdata-new",
)

shot_settings = ShotSettings(
    # logging
    log_settings=LogSettings(
        log_file_path=None,
        file_log_level=logging.WARNING,
        log_file_write_mode="w",
        log_to_console=True,
        console_log_level=logging.WARNING,
        use_custom_logging=False,
    ),
    # data settings
    existing_data_request=None,
    efit_tree_name="analysis",
    # method selection
    run_methods=[],
    run_tags=["all"],
    run_columns=[],
    only_requested_columns=False,
    shot_data_requests=[],
    # timebase settings
    set_times_request="efit",  # use efit timebase
    signal_domain="full",
    use_existing_data_timebase=False,
    interpolation_method="linear",
)

shot_data = handler.get_shots_data(
    shot_ids_request=-1,  # no default value
    shot_settings=shot_settings,
    output_type_request="list",  # output a list of dataframes
    num_processes=1,
)
