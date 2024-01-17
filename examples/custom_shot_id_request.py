from typing import List
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.shot_ids_request import ShotIdsRequest, ShotIdsRequestParams
from disruption_py.settings.shot_settings import ShotSettings

# Create the ShotIdsRequest class that handles getting the shot numbers for the request
# Shot
class CustomShotIdRequest(ShotIdsRequest):
    
    # the _get_shot_ids function takes a set of parameters including a reference to the sql database and returns a list of shot numbers
    def _get_shot_ids(self, params: ShotIdsRequestParams) -> List:
        sql_shot_nums = params.database.query("SELECT shot FROM good_shots WHERE EXTRACT(YEAR FROM entered) BETWEEN 2019 AND 2021;")['shot']
        return sql_shot_nums + [1160405002, 1140523021, 1140523026, 1160620011]


cmod_handler = CModHandler()
shot_settings = ShotSettings(
    set_times_request="efit", # use the efit timebase preset for set_times_request
    run_tags=[],
    run_methods=["_get_ip_parameters"], # only run thr get_ip_parameters method
    output_type_request="ip_data.csv", # automatically uses the CSVOutputRequest preset because of the .csv file descriptor
)
shot_data = cmod_handler.get_shots_data(
    shot_id_request=CustomShotIdRequest(),
    shot_settings=shot_settings,
    num_processes = 4,
)