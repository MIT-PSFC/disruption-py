#!/usr/bin/env python3

"""
The main entrypoint for retrieving DisruptionPy data. 
"""

import logging
from typing import Any, Callable

from disruption_py.core.multiprocess import MultiprocessingShotRetriever
from disruption_py.core.retrieval_manager import RetrievalManager
from disruption_py.core.utils.misc import without_duplicates
from disruption_py.inout.mds import ProcessMDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
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

logger = logging.getLogger("disruption_py")


def get_shots_data(
    shotlist_setting: ShotlistSettingType,
    tokamak: Tokamak = None,
    database_initializer: Callable[..., ShotDatabase] = None,
    mds_connection_initializer: Callable[..., ProcessMDSConnection] = None,
    retrieval_settings: RetrievalSettings = None,
    output_setting: OutputSetting = "list",
    num_processes: int = 1,
    log_settings: LogSettings = None,
) -> Any:
    """
    Get shot data for all shots from shotlist_setting from CMOD.

    Attributes
    ----------
    shotlist_setting : ShotlistSettingType
        Data retrieved for all shotlist specified by the setting. See ShotlistSetting
        for more details.
    retrieval_settings : RetrievalSettings
        The settings that each shot uses when retrieving data. See RetrievalSettings
        for more details. If None, the default values of each setting in
        RetrievalSettings is used.
    output_setting : OutputSetting
        The output type setting to be used when outputting the retrieved data for
        each shot. Note that data is streamed to the output type setting object
        as it is retrieved. Can pass any OutputSettingType that resolves to an
        OutputSetting. See OutputSetting for more details. Defaults to "list".
    num_processes : int
        The number of processes to use for data retrieval. If 1, the data is retrieved
        in serial. If > 1, the data is retrieved in parallel.
    log_settings : LogSettings
        Settings for logging.
    Returns
    -------
    Any
        The value of OutputSetting.get_results. See OutputSetting for more details.
    """
    (log_settings or LogSettings()).setup_logging()

    tokamak = resolve_tokamak_from_environment(tokamak)

    database_initializer = database_initializer or (
        lambda: get_database(tokamak=tokamak)
    )
    database = database_initializer()
    mds_connection_initializer = mds_connection_initializer or (
        lambda: ProcessMDSConnection.from_config(tokamak=tokamak)
    )
    # Clean-up parameters
    if retrieval_settings is None:
        retrieval_settings = RetrievalSettings()

    retrieval_settings.resolve()
    output_setting = resolve_output_setting(output_setting)

    # do not spawn unnecessary processes
    shotlist_setting_params = ShotlistSettingParams(database, tokamak, logger)
    shotlist_list = without_duplicates(
        shotlist_setting_runner(shotlist_setting, shotlist_setting_params)
    )

    num_processes = min(num_processes, len(shotlist_list))

    if num_processes > 1:
        shot_retriever = MultiprocessingShotRetriever(
            database=database,
            num_processes=num_processes,
            output_setting=output_setting,
            retrieval_manager_initializer=(
                lambda: RetrievalManager(
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
            retrieval_settings=retrieval_settings,
            await_complete=True,
        )
    else:
        mds_connection = mds_connection_initializer()
        retrieval_manager = RetrievalManager(
            tokamak=tokamak,
            process_database=database,
            process_mds_conn=mds_connection,
        )
        for shot_id in shotlist_list:
            shot_data = retrieval_manager.get_shot_data(
                shot_id=shot_id,
                retrieval_settings=retrieval_settings,
            )
            if shot_data is None:
                logger.warning(
                    "Not outputting data for shot %s due, data is None.", shot_id
                )
            else:
                output_setting.output_shot(
                    OutputSettingParams(
                        shot_id=shot_id,
                        result=shot_data,
                        database=database,
                        tokamak=tokamak,
                        logger=logger,
                    )
                )

    finish_output_type_setting_params = CompleteOutputSettingParams(
        tokamak=tokamak, logger=logger
    )
    results = output_setting.get_results(finish_output_type_setting_params)
    output_setting.stream_output_cleanup(finish_output_type_setting_params)
    return results


def get_database(
    tokamak: Tokamak = None,
) -> ShotDatabase:
    """
    Get the shot database for the tokamak.
    """
    tokamak = resolve_tokamak_from_environment(tokamak)
    return ShotDatabase.from_config(tokamak=tokamak)


def get_mdsplus_class(
    tokamak: Tokamak = None,
) -> ProcessMDSConnection:
    """
    Get the MDSplus connection for the tokamak.
    """
    tokamak = resolve_tokamak_from_environment(tokamak)
    return ProcessMDSConnection.from_config(tokamak=tokamak)
