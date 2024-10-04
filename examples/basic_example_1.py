#!/usr/bin/env python3

"""
Example usage of `get_shots_data` for retrieving three shots and outputting data 
to csv.
"""

from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data

retrieval_settings = RetrievalSettings(
    # uses the efit timebase when returning data
    time_setting="efit",
    # run all available methods
    run_tags=["all"],
)
shot_data = get_shots_data(
    tokamak="cmod",
    # Retrieve data for the desired shots
    shotlist_setting=[1150805012, 1150805013, 1150805014],
    retrieval_settings=retrieval_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_setting="data.csv",
    num_processes=1,
)
