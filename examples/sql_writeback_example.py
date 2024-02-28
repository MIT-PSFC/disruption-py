from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.settings.output_type_request import SQLOutputRequest

cmod_handler = CModHandler()
shot_settings = ShotSettings(
    # uses the efit timebase when returning data 
    set_times_request="efit",
    efit_tree_name="efit18",
    # run all available methods
    run_tags=["all"],
)
shot_data = cmod_handler.get_shots_data(
    # Retrieve data for the desired shots
    shot_ids_request=[
        1140819002,
		1140819003,
		1140819004,
 	],
    shot_settings=shot_settings,
    
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request=SQLOutputRequest(table_name="disruption_warning_test"),
    
    num_processes = 1
)