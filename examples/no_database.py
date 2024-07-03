#!/usr/bin/env python3

import logging
import os

from disruption_py.io.sql import DummyDatabase
from disruption_py.workflow import get_shots_data
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.retrieval_settings import RetrievalSettings

shot_settings = RetrievalSettings(
    # uses the efit timebase when returning data
    time_setting="ip",
    # run all available methods
    run_tags=["all"],
)
shot_data = get_shots_data(
    tokamak="d3d",
    # Retrieve data for the desired shots
    shotlist_setting=[161228, 161237, 166177, 166253],
    shot_settings=shot_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_setting=f"/tmp/{os.environ['USER']}/data.csv",
    num_processes=1,
    database_initializer=DummyDatabase.initializer,  # use dummy database
    log_settings=LogSettings(
        console_log_level=logging.DEBUG,
    ),
)
