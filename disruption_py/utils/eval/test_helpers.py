from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings import LogSettings, ShotSettings


import pandas as pd


import logging
from typing import Dict, List

from disruption_py.utils.constants import TIME_CONST


def get_mdsplus_data(cmod_handler : CModHandler, shot_ids : List[int]) -> Dict[int, pd.DataFrame]:
    """
    Get MDSplus data for a list of shots.

     Returns
    -------
    Dict[int, pd.DataFrame]
        Dictionary mapping shot IDs to retrieved MDSplus data.
      """
    shot_settings = ShotSettings(
        efit_tree_name="efit18",
        set_times_request="efit",
        log_settings=LogSettings(
            log_to_console=False,
            log_file_path="tests/cmod.log",
            log_file_write_mode="w",
            file_log_level=logging.DEBUG
        )
    )
    shot_data = cmod_handler.get_shots_data(
        shot_ids_request=shot_ids,
        shot_settings=shot_settings,
        output_type_request="dict",
    )
    return shot_data


def get_sql_data_for_mdsplus(cmod_handler : CModHandler, shot_ids : List[int], mdsplus_data : Dict[int, pd.DataFrame]) -> Dict[int, pd.DataFrame]:
    """
    Get SQL data for a list of shots and map onto the timebase of the supplied MDSplus data.

    Returns
    -------
    Dict[int, pd.DataFrame]
        Dictionary mapping shot IDs to retrieved SQL data.
    """
    shot_data = {}
    for shot_id in shot_ids:
        times = mdsplus_data[shot_id]['time']
        sql_data = cmod_handler.database.get_shots_data([shot_id])
        shot_data[shot_id] = pd.merge_asof(times.to_frame(), sql_data, on='time', direction='nearest', tolerance=TIME_CONST)
    return shot_data