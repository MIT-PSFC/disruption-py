from disruption_py.main import get_database, get_shots_data
from disruption_py.settings.output_type_request import SQLOutputRequest
from disruption_py.settings.shot_settings import ShotSettings

shot_settings = ShotSettings(
    # uses the efit timebase when returning data
    set_times_request="efit",
    efit_tree_name="efit18",
    # run all available methods
    run_tags=["all"],
)
shot_ids = [1140819005, 1140819009]
shot_data = get_shots_data(
    # Retrieve data for the desired shots
    shot_ids_request=shot_ids,
    shot_settings=shot_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request=SQLOutputRequest(table_name="disruption_warning_test"),
    num_processes=1,
)


cmod_database = get_database(tokamak="cmod")
result = cmod_database.get_shots_data(shot_ids, sql_table="disruption_warning_test")
print(result)
