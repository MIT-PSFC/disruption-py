import logging
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.set_times_request import SetTimesRequest, SetTimesRequestParams
from disruption_py.settings.shot_data_request import ShotDataRequest, ShotDataRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.method_caching import parameter_cached_method


cmod_handler = CModHandler()

class PRadTime(SetTimesRequest): 
    def _get_times(self, params : SetTimesRequestParams):
        time_array, = params.mds_conn.get_dims(r"\twopi_diode", tree_name='spectroscopy')
        time_array = time_array[time_array > 0]
        return time_array

shot_settings = ShotSettings(
    set_times_request=PRadTime(),    
    run_tags=[],
	run_columns=["p_rad"],
    only_requested_columns=True,
	signal_domain="flattop",
	log_settings=LogSettings(
		log_to_console=False,
		log_file_path="examples/last_log.log",
	)
)

shot_data = cmod_handler.get_shots_data(
    # use the shot ids listed in the shot_ids_of_interest.txt file
    shot_ids_request="cmod_non_disruptions_ids_not_blacklist_mini",
    shot_settings=shot_settings,
    
    # automatically stream retrieved data to a csv file by passing in a file path ending in .csv
    # this works by automatically using the CSVOutputRequest preset with the passed filename
    output_type_request="examples/p_rad.csv", 
    
    # retrieve data quickly using 4 processes
    num_processes = 4,
)