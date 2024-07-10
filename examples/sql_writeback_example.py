#!/usr/bin/env python3

from disruption_py.workflow import get_database, get_shots_data
from disruption_py.settings.output_setting import SQLOutputSetting
from disruption_py.settings.retrieval_settings import RetrievalSettings

retrieval_settings = RetrievalSettings(
    # uses the efit timebase when returning data
    time_setting="efit",
    efit_tree_name="efit18",
    # run all available methods
    run_tags=["all"],
)
shotlist = [1140819005, 1140819009]
shot_data = get_shots_data(
    # Retrieve data for the desired shots
    shotlist_setting=shotlist,
    retrieval_settings=retrieval_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_setting=SQLOutputSetting(table_name="disruption_warning_test"),
    num_processes=1,
)


cmod_database = get_database(tokamak="cmod")
result = cmod_database.get_shots_data(shotlist, sql_table="disruption_warning_test")
print(result)
