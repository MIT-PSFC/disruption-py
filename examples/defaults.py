#!/usr/bin/env python3

"""
Example usage of `get_shots_data` with all the default arguments explicitly assigned.
"""

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

retrieval_settings = RetrievalSettings(
    efit_nickname_setting="disruption",
    # method/column selection
    # default None: all methods/columns
    run_methods=None,
    run_columns=None,
    only_requested_columns=False,
    custom_physics_methods=[],
    # timebase settings
    time_setting="disruption_warning",
    domain_setting="full",
)

shot_data = get_shots_data(
    # required argument
    shotlist_setting=[],
    # default None: detect from environment
    tokamak=None,
    # default None: standard SQL/MDSplus connection
    database_initializer=None,
    mds_connection_initializer=None,
    retrieval_settings=retrieval_settings,
    output_setting="dataset",
    num_processes=1,
    log_settings=LogSettings(
        # default None: "output.log" in temporary session folder
        file_path=None,
        file_level="DEBUG",
        # default None: VERBOSE, or higher based on number of shots
        console_level=None,
    ),
)
