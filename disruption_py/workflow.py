#!/usr/bin/env python3

import logging
from typing import Any, Callable

from disruption_py.io.sql import ShotDatabase
from disruption_py.io.mds import ProcessMDSConnection
from disruption_py.settings import Settings
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.output_setting import (
    CompleteOutputSettingParams,
    OutputSetting,
    OutputSettingParams,
    resolve_output_setting,
)
from disruption_py.settings.shotlist_setting import (
    ShotlistSettingParams,
    ShotlistSettingType,
    shotlist_setting_runner,
)
from disruption_py.machine.tokamak import Tokamak, resolve_tokamak
from disruption_py.machine.factory import (
    get_tokamak_shot_manager,
)
from disruption_py.utils.multiprocessing_helper import MultiprocessingShotRetriever
from disruption_py.utils.utils import without_duplicates

logger = logging.getLogger("disruption_py")


def get_shots_data(
    shotlist_request: ShotlistSettingType,
    tokamak: Tokamak = None,
    database_initializer: Callable[..., ShotDatabase] = None,
    mds_connection_initializer: Callable[..., ProcessMDSConnection] = None,
    shot_settings: Settings = None,
    output_type_request: OutputSetting = "list",
    num_processes: int = 1,
    log_settings: LogSettings = None,
) -> Any:
    """
    Get shot data for all shots from shotlist_request from CMOD.

    Attributes
    ----------
    shotlist_request : ShotIdsRequestType
        Data retrieved for all shotlist specified by the request. See ShotIdsRequest for more details.
    shot_settings : ShotSettings
        The settings that each shot uses when retrieving data. See ShotSettings for more details.
        If None, the default values of each setting in ShotSettings is used.
    output_type_request : OutputTypeRequest
        The output type request to be used when outputting the retrieved data for each shot. Note that data
        is streamed to the output type request object as it is retrieved. Can pass any OutputTypeRequestType
        that resolves to an OutputTypeRequest. See OutputTypeRequest for more details. Defaults to "list".
    num_processes : int
        The number of processes to use for data retrieval. If 1, the data is retrieved in serial.
        If > 1, the data is retrieved in parallel.
    log_settings : LogSettings
        Settings for logging.
    Returns
    -------
    Any
        The value of OutputTypeRequest.get_results, where OutputTypeRequest is specified in
        shot_settings. See OutputTypeRequest for more details.
    """
    (log_settings or LogSettings.default()).setup_logging()

    tokamak = resolve_tokamak(tokamak)

    database_initializer = database_initializer or (
        lambda: get_database(tokamak=tokamak)
    )
    database = database_initializer()
    mds_connection_initializer = mds_connection_initializer or (
        lambda: ProcessMDSConnection.from_config(tokamak=tokamak)
    )
    # Clean-up parameters
    if shot_settings is None:
        shot_settings = Settings()

    shot_settings.resolve()
    output_type_request = resolve_output_setting(output_type_request)

    # do not spawn unnecessary processes
    shot_manager_cls = get_tokamak_shot_manager(tokamak)

    shotlist_request_params = ShotlistSettingParams(database, tokamak, logger)
    shotlist_list = without_duplicates(
        shotlist_setting_runner(shotlist_request, shotlist_request_params)
    )

    num_processes = min(num_processes, len(shotlist_list))

    if num_processes > 1:
        shot_retriever = MultiprocessingShotRetriever(
            database=database,
            num_processes=num_processes,
            output_type_request=output_type_request,
            shot_manager_initializer=(
                lambda: shot_manager_cls(
                    tokamak=tokamak,
                    process_database=database_initializer(),
                    process_mds_conn=mds_connection_initializer(),
                )
            ),
            tokamak=tokamak,
            logger=logger,
        )
        shot_retriever.run(
            shotlist_list=shotlist_list,
            shot_settings=shot_settings,
            await_complete=True,
        )
    else:
        mds_connection = mds_connection_initializer()
        shot_manager = shot_manager_cls(
            tokamak=tokamak,
            process_database=database,
            process_mds_conn=mds_connection,
        )
        for shot_id in shotlist_list:
            shot_data = shot_manager.get_shot_data(
                shot_id=shot_id,
                shot_settings=shot_settings,
            )
            if shot_data is None:
                logger.warning(
                    f"Not outputting data for shot {shot_id} due, data is None."
                )
            else:
                output_type_request.output_shot(
                    OutputSettingParams(
                        shot_id=shot_id,
                        result=shot_data,
                        database=database,
                        tokamak=tokamak,
                        logger=logger,
                    )
                )

    finish_output_type_request_params = CompleteOutputSettingParams(
        tokamak=tokamak, logger=logger
    )
    results = output_type_request.get_results(finish_output_type_request_params)
    output_type_request.stream_output_cleanup(finish_output_type_request_params)
    return results


def get_database(
    tokamak: Tokamak = None,
) -> ShotDatabase:
    """
    Get the shot database for the tokamak.
    """
    tokamak = resolve_tokamak(tokamak)
    return ShotDatabase.from_config(tokamak=tokamak)


def get_mdsplus_class(
    tokamak: Tokamak = None,
) -> ProcessMDSConnection:
    """
    Get the MDSplus connection for the tokamak.
    """
    tokamak = resolve_tokamak(tokamak)
    return ProcessMDSConnection.from_config(tokamak=tokamak)
