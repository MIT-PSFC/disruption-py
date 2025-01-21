#!/usr/bin/env python3

"""
Example usage of `get_shots_data` demonstrating how to provide a customized 
shotlist_setting by subclassing `ShotlistSetting`.
"""

from typing import List

from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.settings.shotlist_setting import (
    ShotlistSetting,
    ShotlistSettingParams,
)
from disruption_py.workflow import get_shots_data


class CustomShotlistSetting(ShotlistSetting):
    """The ShotlistSetting class handles getting the needed shotlist"""

    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        """
        Takes a set of parameters including a reference to the sql database and
        return a list of shot numbers.
        """
        sql_shot_nums = params.database.query(
            "SELECT shot FROM good_shots WHERE EXTRACT(YEAR FROM entered) BETWEEN 2019 AND 2021;"
        )["shot"]
        return sql_shot_nums + [1160405002, 1140523021, 1140523026, 1160620011]


retrieval_settings = RetrievalSettings(
    # use the efit timebase preset for time_setting
    time_setting="efit",
    # only run thr get_ip_parameters method
    run_methods=["_get_ip_parameters"],
)
shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=CustomShotlistSetting(),
    retrieval_settings=retrieval_settings,
    # automatically uses the CSVOutputSetting preset because of the .csv file descriptor
    output_setting="ip_data.csv",
    num_processes=4,
)
