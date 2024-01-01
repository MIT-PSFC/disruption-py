from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings import ShotSettings

cmod_handler = CModHandler()
shot_settings = ShotSettings(
	existing_data_request="sql",
	set_times_request="magnetics004",
	run_tags=[],
	run_methods=["_get_ip_parameters"],
	output_type_request="ip_data.csv",
)
cmod_handler.get_shots_data(
	shot_ids_request=[1160405002, 1140523021, 1140523026, 1160620011],
	shot_settings=shot_settings,
	num_processes = 4,
)