from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.shot_settings import ShotSettings


cmod_handler = CModHandler()
shot_settings = ShotSettings(
    # retrieve existing data from the disruption_warnings table use it to prepopulate queries
    existing_data_request="sql",
    
    # use the efit timebase preset for set_times_request
    # uses the efit timebase when returning data 
    set_times_request="efit",
        
    # only run the get_ip_parameters method
    run_methods=["_get_ip_parameters"],
    run_tags=[],
    
    # automatically uses the CSVOutputRequest preset because of the .csv file descriptor
    # streams outputted data to the ip_data.csv file
    output_type_request="ip_data.csv", 
)
shot_data = cmod_handler.get_shots_data(
    shot_ids_request="shot_ids_of_interest.txt", # use the shot ids listed in the 
    shot_settings=shot_settings,
    num_processes = 4,
)