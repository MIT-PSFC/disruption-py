#!/usr/bin/env python3

from typing import List

from disruption_py.workflow import get_shots_data
from disruption_py.settings.shotlist_setting import (
    ShotlistSetting,
    ShotlistSettingParams,
)
from disruption_py.settings.settings import Settings


# Create the ShotIdsRequest class that handles getting the shot numbers for the request
# Shot
class CustomShotIdRequest(ShotlistSetting):

    # the _get_shotlist function takes a set of parameters including a reference to the sql database and returns a list of shot numbers
    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        sql_shot_nums = params.database.query(
            "SELECT shot FROM good_shots WHERE EXTRACT(YEAR FROM entered) BETWEEN 2019 AND 2021;"
        )["shot"]
        return sql_shot_nums + [1160405002, 1140523021, 1140523026, 1160620011]


shot_settings = Settings(
    # use the efit timebase preset for set_times_request
    set_times_request="efit",
    run_tags=[],
    # only run thr get_ip_parameters method
    run_methods=["_get_ip_parameters"],
)
shot_data = get_shots_data(
    tokamak="cmod",
    shot_id_request=CustomShotIdRequest(),
    shot_settings=shot_settings,
    # automatically uses the CSVOutputRequest preset because of the .csv file descriptor
    output_type_request="ip_data.csv",
    num_processes=4,
)
