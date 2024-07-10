#!/usr/bin/env python3

import logging

from disruption_py.workflow import get_shots_data
from disruption_py.settings import LogSettings, RetrievalSettings

retrieval_settings = RetrievalSettings(
    # data settings
    input_setting=None,
    efit_tree_name="analysis",
    # method selection
    run_methods=[],
    run_tags=["all"],
    run_columns=[],
    only_requested_columns=False,
    custom_physics_methods=[],
    # timebase settings
    time_setting="disruption_warning",  # use efit timebase
    domain_setting="full",
    use_input_setting_timebase=False,
    interpolation_method="linear",
)

shot_data = get_shots_data(
    tokamak=None,  # defaults to tokamak value detected from environement
    shotlist_setting=-1,  # no default value
    database_initializer=None,  # defaults to connection for tokamak
    mds_connection_initializer=None,  # defaults to mds plus server string for tokamak
    retrieval_settings=retrieval_settings,
    output_setting="list",  # output a list of dataframes
    num_processes=1,
    log_settings=LogSettings(  # logging
        log_file_path=None,
        file_log_level=logging.WARNING,
        log_file_write_mode="w",
        log_to_console=True,
        console_log_level=logging.WARNING,
        use_custom_logging=False,
    ),
)
