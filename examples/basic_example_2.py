#!/usr/bin/env python3

from disruption_py.workflow import get_shots_data
from disruption_py.settings.retrieval_settings import RetrievalSettings

shot_settings = RetrievalSettings(
    # retrieve existing data from the disruption_warnings table use it to prepopulate queries
    input_setting="sql",
    # use the efit timebase preset for the time_setting
    # uses the efit timebase when returning data
    time_setting="efit",
    # run the get_ip_parameters method
    run_methods=["_get_ip_parameters"],
    run_tags=[],
    # run the method that  returns the ne_peaking column
    run_columns=["ne_peaking"],
)
shot_data = get_shots_data(
    tokamak="cmod",
    # use the shotlist listed in the shotlist_of_interest.txt file
    shotlist_setting="shotlist_of_interest.txt",
    shot_settings=shot_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    # this works by automatically using the CSVOutputSetting preset with the passed filename
    output_setting="ip_data.csv",
    # retrieve data quickly using 4 processes
    num_processes=4,
)
