#!/usr/bin/env python3

"""
Example usage of `get_shots_data` with all the default arguments explicitly assigned.
"""

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

retrieval_settings = RetrievalSettings(
    # data settings
    cache_setting=None,
    efit_nickname_setting="disruption",  # defaults to use the "disruption" efit tree for disruptive shots
    # method selection.
    run_methods=None,  # defaults to run all available methods for the machine
    run_columns=None,  # defaults to return all available columns for the machine
    only_requested_columns=False,  # defaults to return all columns retrieved by the methods
    custom_physics_methods=[],
    # timebase settings
    time_setting="disruption_warning",  # defaults to return data in the "disruption_warning" time base
    domain_setting="full",  # defaults to getting the entire shot
    use_cache_setting_timebase=False,
    interpolation_method="linear",  # defaults to use linear interpolation
)

shot_data = get_shots_data(
    shotlist_setting=[],  # required
    tokamak=None,  # defaults to detect from environment
    database_initializer=None,  # defaults to SQL connection for tokamak
    mds_connection_initializer=None,  # defaults to MDSplus server string for tokamak
    retrieval_settings=retrieval_settings,
    output_setting="dataframe",
    num_processes=1,
    log_settings=LogSettings(
        log_file_path=None,  # defaults to "output.log" in tmp folder for the session
        file_log_level="DEBUG",
        log_file_write_mode="w",
        log_to_console=True,
        console_log_level=None,  # defaults to VERBOSE but varies based on number of shots
        use_custom_logging=False,
    ),
)
