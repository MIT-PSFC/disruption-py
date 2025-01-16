#!/usr/bin/env python3

"""
Example usage of `get_shots_data` with all the default arguments explicitly assigned.
"""

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

retrieval_settings = RetrievalSettings(
    # data settings
    cache_setting=None,
    efit_nickname_setting="disruption",
    # method selection
    run_methods=[],
    run_tags=["all"],
    run_columns=[],
    only_requested_columns=False,
    custom_physics_methods=[],
    # timebase settings
    time_setting="disruption_warning",
    domain_setting="full",
    use_cache_setting_timebase=False,
    interpolation_method="linear",
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
