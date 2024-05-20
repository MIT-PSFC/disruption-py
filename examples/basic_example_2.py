from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.shot_settings import ShotSettings


cmod_handler = CModHandler()
shot_settings = ShotSettings(
    # retrieve existing data from the disruption_warnings table use it to prepopulate queries
    existing_data_request="sql",
    # use the efit timebase preset for set_times_request
    # uses the efit timebase when returning data
    set_times_request="efit",
    # run the get_ip_parameters method
    run_methods=["_get_ip_parameters"],
    run_tags=[],
    # run the method that  returns the ne_peaking column
    run_columns=["ne_peaking"],
)
shot_data = cmod_handler.get_shots_data(
    # use the shot ids listed in the shot_ids_of_interest.txt file
    shot_ids_request="shot_ids_of_interest.txt",
    shot_settings=shot_settings,
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    # this works by automatically using the CSVOutputRequest preset with the passed filename
    output_type_request="ip_data.csv",
    # retrieve data quickly using 4 processes
    num_processes=4,
)
