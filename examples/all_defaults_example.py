#!/usr/bin/env python3

import logging

from disruption_py.main import get_shots_data
from disruption_py.settings import LogSettings, ShotSettings

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
    custom_parameter_methods=[],
    # timebase settings
    set_times_request="disruption_warning",  # use efit timebase
    signal_domain="full",
    use_existing_data_timebase=False,
    interpolation_method="linear",
)

shot_data = get_shots_data(
    tokamak=None,  # defaults to tokamak value detected from environement
    shot_ids_request=-1,  # no default value
    database_initializer=None,  # defaults to connection for tokamak
    mds_connection_str=None,  # defaults to mds plus server string for tokamak
    shot_settings=shot_settings,
    output_type_request="list",  # output a list of dataframes
    num_processes=1,
)
