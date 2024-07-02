#!/usr/bin/env python3

from disruption_py.workflow import get_shots_data
from disruption_py.settings.shot_settings import ShotSettings

shot_settings = ShotSettings(
    # uses the efit timebase when returning data
    set_times_request="efit",
    # run all available methods
    run_tags=["all"],
)
shot_data = get_shots_data(
    tokamak="cmod",
    # Retrieve data for the desired shots
    shot_ids_request=[1150805012, 1150805013, 1150805014],
    shot_settings=shot_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request="data.csv",
    num_processes=1,
)
