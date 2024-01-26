import logging
import sys

sys.path.append("/home/joshlor/Documents/disruption-py/")
from disruption_py.settings.log_settings import LogSettings
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.shot_settings import ShotSettings

cmod_handler = CModHandler()
shot_settings = ShotSettings(
    # uses the efit timebase when returning data 
    set_times_request="efit",
    
    # run all available methods
    run_tags=["all"],
    
    efit_tree_name="efit18",
    
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    
)

shot_data = cmod_handler.get_shots_data(
    shot_ids_request=[1160405002, 1140523021, 1140523026, 1160620011], # Retrieve data for the desired shots
    shot_settings=shot_settings,
    
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request="test/ip_data.hdf5", 
    
    num_processes=1,
)

