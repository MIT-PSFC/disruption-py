import sys
sys.path.append("/home/lorinczj/disruption-py")

import logging
import os
from disruption_py.settings.log_settings import LogSettings
from disruption_py.handlers.d3d_handler import D3DHandler
from disruption_py.databases.dummy_database import DummyDatabase
from disruption_py.settings.shot_settings import ShotSettings

USER = os.getenv('USER')

cmod_handler = D3DHandler(database_initializer=DummyDatabase.default)
shot_settings = ShotSettings(
    # uses the efit timebase when returning data 
    set_times_request="ip",
    
    # run all available methods
    run_tags=["all"],
    
    log_settings=LogSettings(
        console_log_level=logging.DEBUG,
	),
)
shot_data = cmod_handler.get_shots_data(
    # Retrieve data for the desired shots
    shot_ids_request=[161228, 161237, 166177, 166253],
    shot_settings=shot_settings,
    
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    output_type_request=f"/cscratch/{USER}/data.csv",
    
    num_processes = 1
)