#!/usr/bin/env python3

import logging
import os

from disruption_py.database import DummyDatabase
from disruption_py.main import get_shots_data
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.shot_settings import ShotSettings

shot_settings = ShotSettings(
    # uses the efit timebase when returning data
    set_times_request="ip",
    # run all available methods
    run_tags=["all"],
    log_settings=LogSettings(
        console_log_level=logging.DEBUG,
    ),
)
shot_data = get_shots_data(
    tokamak="d3d",
    # Retrieve data for the desired shots
    shot_ids_request=[161228, 161237, 166177, 166253],
    shot_settings=shot_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request=f"/tmp/{os.environ['USER']}/data.csv",
    num_processes=1,
    database_initializer=DummyDatabase.initializer,  # use dummy database
)
