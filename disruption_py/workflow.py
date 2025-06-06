#!/usr/bin/env python3

"""
main workflow
"""

import argparse
import time
from itertools import repeat
from multiprocessing import Pool
from typing import Any, Callable

from loguru import logger
from tqdm.auto import tqdm

from disruption_py.config import config
from disruption_py.core.retrieval_manager import RetrievalManager
from disruption_py.core.utils.misc import (
    get_elapsed_time,
    without_duplicates,
)
from disruption_py.inout.mds import ProcessMDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.settings.log_settings import LogSettings, resolve_log_settings
from disruption_py.settings.output_setting import (
    OutputSetting,
    OutputSettingParams,
    resolve_output_setting,
)
from disruption_py.settings.shotlist_setting import (
    ShotlistSettingParams,
    ShotlistSettingType,
    shotlist_setting_runner,
)


def _execute_retrieval(args):
    """
    Wrapper around getting shot data for a single shot to ensure the arguments
    are all multiprocessing compatible (e.g. no lambdas passed as args).

    Params
    ------
    args : List
        tokamak, database initializer, mds connection initializer, retrieval
        settings, and the shot id

    Returns
    -------
    tuple of shot id and the dataframe
    """
    tokamak, db_init, mds_init, retrieval_settings, shot_id = args
    database = _get_database_instance(tokamak, db_init)
    mds_conn = _get_mds_instance(tokamak, mds_init)

    retrieval_manager = RetrievalManager(
        tokamak=tokamak,
        process_database=database,
        process_mds_conn=mds_conn,
    )
    return shot_id, retrieval_manager.get_shot_data(shot_id, retrieval_settings)


def get_shots_data(
    shotlist_setting: ShotlistSettingType,
    tokamak: Tokamak = None,
    database_initializer: Callable[..., ShotDatabase] = None,
    mds_connection_initializer: Callable[..., ProcessMDSConnection] = None,
    retrieval_settings: RetrievalSettings = None,
    output_setting: OutputSetting = "dataset",
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
        OutputSetting. See OutputSetting for more details. Defaults to "dataset".
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
    log_settings = resolve_log_settings(log_settings)
    log_settings.setup_logging()

    tokamak = resolve_tokamak_from_environment(tokamak)
    database = _get_database_instance(tokamak, database_initializer)
    # Clean-up parameters
    if retrieval_settings is None:
        retrieval_settings = RetrievalSettings()

    retrieval_settings.resolve()
    output_setting = resolve_output_setting(output_setting)

    # do not spawn unnecessary processes
    shotlist_setting_params = ShotlistSettingParams(database, tokamak)
    shotlist_list = without_duplicates(
        shotlist_setting_runner(shotlist_setting, shotlist_setting_params)
    )
    num_processes = min(num_processes, len(shotlist_list))
    if num_processes < 1:
        logger.critical("Nothing to do!")
        return None

    # Dynamically set the console log level based on the number of shots
    if log_settings.console_level is None:
        log_settings.reset_handlers(num_shots=len(shotlist_list))

    # log start
    logger.info(
        "Starting workflow: {n:,} shot{s} / {m} process{p}",
        n=len(shotlist_list),
        s="s" if len(shotlist_list) > 1 else "",
        m=num_processes,
        p="es" if num_processes > 1 else "",
    )

    took = -time.time()
    with Pool(processes=num_processes) as pool:
        args = zip(
            repeat(tokamak),
            repeat(database_initializer),
            repeat(mds_connection_initializer),
            repeat(retrieval_settings),
            shotlist_list,
        )
        num_success = 0

        for shot_id, shot_data in tqdm(
            pool.imap(_execute_retrieval, args), total=len(shotlist_list), leave=False
        ):
            if shot_data is not None:
                num_success += 1
                output_setting.output_shot(
                    OutputSettingParams(
                        shot_id=shot_id,
                        result=shot_data,
                        tokamak=tokamak,
                    )
                )
    took += time.time()

    # log stop
    total = len(shotlist_list)
    percent_success = num_success / total * 100
    logger.info(
        "Completed workflow: "
        "retrieved {num_success:,}/{total:,} shots ({percent_success:.2f}%) "
        "in {elapsed} ({each:.3f} s/shot)",
        num_success=num_success,
        total=total,
        percent_success=percent_success,
        elapsed=get_elapsed_time(took),
        each=took / total,
    )

    results = output_setting.get_results()
    output_setting.to_disk()
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


def _get_database_instance(tokamak, database_initializer):
    """
    Create database instance
    """
    if database_initializer:
        return database_initializer()
    return get_database(tokamak)


def _get_mds_instance(tokamak, mds_connection_initializer):
    """
    Create MDSplus instance
    """
    if mds_connection_initializer:
        return mds_connection_initializer()
    return get_mdsplus_class(tokamak)


def run(tokamak, methods, shots, efit_tree, time_base, output, processes, log_level):
    """
    simple workflow.
    """

    if not tokamak:
        tokamak = resolve_tokamak_from_environment()
    if not shots:
        shots, *_ = config(tokamak).tests.shots.values()
    methods = methods or None

    sett = RetrievalSettings(
        efit_nickname_setting=efit_tree, time_setting=time_base, run_methods=methods
    )

    return get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shots,
        retrieval_settings=sett,
        num_processes=processes,
        log_settings=log_level,
        output_setting=output,
    )


def cli():
    """
    simple argument parser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("shots", nargs="*")
    parser.add_argument("-t", "--tokamak", type=str)
    parser.add_argument("-m", "--methods", type=str, action="append")
    parser.add_argument("-e", "--efit-tree", type=str, default="disruption")
    parser.add_argument("-b", "--time-base", type=str, default="disruption_warning")
    parser.add_argument("-o", "--output", type=str, default="dataset")
    parser.add_argument("-p", "--processes", type=int, default=1)
    parser.add_argument("-l", "--log-level", type=str, default="VERBOSE")

    return run(**vars(parser.parse_args()))


if __name__ == "__main__":
    out = cli()
    print(out)
