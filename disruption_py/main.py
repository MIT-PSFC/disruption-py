
#!/usr/bin/env python3

import logging
from typing import Any, Callable

from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_shot_manager, get_database_initializer_for_tokamak, get_mds_connection_str_for_tokamak
from disruption_py.utils.multiprocessing_helper import MultiprocessingShotRetriever
from disruption_py.mdsplus_integration.mds_connection import ProcessMDSConnection
from disruption_py.settings import ShotSettings
from disruption_py.settings.output_type_request import (
    FinishOutputTypeRequestParams,
    OutputTypeRequest,
    ResultOutputTypeRequestParams,
    resolve_output_type_request,
)
from disruption_py.settings.shot_ids_request import (
    ShotIdsRequestParams,
    ShotIdsRequestType,
    shot_ids_request_runner,
)
from disruption_py.shots.shot_manager import ShotManager
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.utils import without_duplicates

logger = logging.getLogger("disruption_py")

class ConnectionHandler():
    """
    Class used to handle the connection to the MDSplus server.
    """
    def __init__(
        self,
        tokamak: Tokamak,
        database_initializer: Callable[..., ShotDatabase],
        mds_connection_str,
    ):
        self.database_initializer = get_database_initializer_for_tokamak(tokamak, database_initializer)
        self.mds_connection_initializer = lambda: ProcessMDSConnection(
            get_mds_connection_str_for_tokamak(tokamak, mds_connection_str)
        )

    @property
    def database(self) -> ShotDatabase:
        """Reference to the sql shot logbook.

        Returns
        -------
        D3DDatabase
            Database object with an open connection to the D3D sql database.
        """
        if not hasattr(self, "_database"):
            self._database = self.database_initializer()
        return self._database

    @property
    def mds_connection(self) -> ProcessMDSConnection:
        """Reference to the MDSplus connection.

        Returns
        -------
        ProcessMDSConnection
            MDSplus connection object for D3D.
        """
        if not hasattr(self, "_mds_connection"):
            self._mds_connection = self.mds_connection_initializer()
        return self._mds_connection


def get_shots_data(
    tokamak: Tokamak,
    shot_ids_request: ShotIdsRequestType,
    database_initializer: Callable[..., ShotDatabase] = None,
    mds_connection_str : str = None,
    shot_settings: ShotSettings = None,
    output_type_request: OutputTypeRequest = "list",
    num_processes: int = 1,
) -> Any:
    """
    Get shot data for all shots from shot_ids_request from CMOD.

    Attributes
    ----------
    shot_ids_request : ShotIdsRequestType
        Data retrieved for all shot_ids specified by the request. See ShotIdsRequest for more details.
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

    Returns
    -------
    Any
        The value of OutputTypeRequest.get_results, where OutputTypeRequest is specified in
        shot_settings. See OutputTypeRequest for more details.
    """

    # Clean-up parameters
    if shot_settings is None:
        shot_settings = ShotSettings()
  
    shot_settings.resolve()
    output_type_request = resolve_output_type_request(output_type_request)
    connection_manager = ConnectionHandler(tokamak, database_initializer, mds_connection_str)

    # do not spawn unnecessary processes
    num_processes = min(num_processes, len(shot_ids_list))
    shot_manager_cls = get_tokamak_shot_manager(tokamak)
    
    shot_ids_request_params = ShotIdsRequestParams(connection_manager.database, tokamak, logger)
    shot_ids_list = without_duplicates(
        shot_ids_request_runner(shot_ids_request, shot_ids_request_params)
    )
    
    if num_processes > 1:
        shot_retriever = MultiprocessingShotRetriever(
            database=connection_manager.database,
            num_processes=num_processes,
            shot_settings=shot_settings,
            output_type_request=output_type_request,
            process_prop_initializers={
                # initialize connections for individual processes
                "shot_manager": (
                    lambda: shot_manager_cls(
                        process_database=connection_manager.database_initializer(),
                        process_mds_conn=connection_manager.mds_connection_initializer(),
                    )
                )
            },
            tokamak=tokamak,
            logger=logger,
        )
        shot_retriever.run(
            shot_creator_f=shot_manager_cls.get_shot_data,
            shot_ids_list=shot_ids_list,
            await_complete=True,
        )
    else:
        shot_manager = shot_manager_cls(
            process_database=connection_manager.database,
            process_mds_conn=connection_manager.mds_connection,
        )
        for shot_id in shot_ids_list:
            shot_data = shot_manager.get_shot_data(
                shot_id=shot_id,
                shot_manager=shot_manager,
                shot_settings=shot_settings,
            )
            if shot_data is None:
                logger.warning(
                    f"Not outputting data for shot {shot_id} due, data is None."
                )
            else:
                output_type_request.output_shot(
                    ResultOutputTypeRequestParams(
                        shot_id=shot_id,
                        result=shot_data,
                        database=connection_manager.database,
                        tokamak=tokamak,
                        logger=logger,
                    )
                )

    finish_output_type_request_params = FinishOutputTypeRequestParams(
        tokamak=tokamak, logger=logger
    )
    results = output_type_request.get_results(finish_output_type_request_params)
    output_type_request.stream_output_cleanup(finish_output_type_request_params)
    return results
