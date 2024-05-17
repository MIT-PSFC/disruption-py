from abc import ABC, abstractmethod
from typing import Callable, Any
import traceback
from disruption_py.databases.database import ShotDatabase
from disruption_py.handlers.multiprocessing_helper import MultiprocessingShotRetriever
from disruption_py.mdsplus_integration.mds_connection import ProcessMDSConnection
from disruption_py.settings.shot_ids_request import (
    ShotIdsRequestParams,
    ShotIdsRequestType,
    shot_ids_request_runner,
)
from disruption_py.settings.output_type_request import (
    OutputTypeRequest,
    ResultOutputTypeRequestParams,
    FinishOutputTypeRequestParams,
    resolve_output_type_request,
)
from disruption_py.settings import ShotSettings
from disruption_py.shots.shot_manager import ShotManager
from disruption_py.utils.mappings.tokamak import Tokamak
import pandas as pd
import logging

from disruption_py.utils.utils import without_duplicates


class Handler(ABC):
    """Base class for all handlers."""

    logger = logging.getLogger("disruption_py")

    def __init__(
        self,
        database_initializer: Callable[..., ShotDatabase],
        mds_connection_str,
        **kwargs,
    ):
        self.database_initializer = database_initializer
        self.mds_connection_initializer = lambda: ProcessMDSConnection(
            mds_connection_str
        )
        self.kwargs = kwargs

    @abstractmethod
    def get_shot_manager_cls(self):
        pass

    @abstractmethod
    def get_tokamak(self):
        pass

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

    @classmethod
    def _get_shot_data(
        cls, shot_id, shot_manager: ShotManager, shot_settings: ShotSettings
    ) -> pd.DataFrame:
        """
        Get data for a single shot. May be run across different processes.
        """
        cls.logger.info(f"starting {shot_id}")
        try:
            shot_props = shot_manager.shot_setup(
                shot_id=int(shot_id),
                shot_settings=shot_settings,
            )
            retrieved_data = shot_manager.shot_data_retrieval(
                shot_props=shot_props, shot_settings=shot_settings
            )
            shot_manager.shot_cleanup(shot_props)
            cls.logger.info(f"completed {shot_id}")
            return retrieved_data
        except Exception as e:
            cls.logger.warning(
                f"[Shot {shot_id}]: fatal error {traceback.format_exc()}"
            )
            cls.logger.error(f"failed {shot_id} with error {e}")
            return None

    def get_shots_data(
        self,
        shot_ids_request: ShotIdsRequestType,
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

        shot_ids_request_params = ShotIdsRequestParams(
            self.database, Tokamak.D3D, self.logger
        )
        shot_ids_list = without_duplicates(
            shot_ids_request_runner(shot_ids_request, shot_ids_request_params)
        )

        # do not spawn unnecessary processes
        num_processes = min(num_processes, len(shot_ids_list))

        if num_processes > 1:
            shot_retriever = MultiprocessingShotRetriever(
                database=self.database,
                num_processes=num_processes,
                shot_settings=shot_settings,
                output_type_request=output_type_request,
                process_prop_initializers={
                    # initialize connections for individual processes
                    "shot_manager": (
                        lambda: self.get_shot_manager_cls()(
                            process_database=self.database_initializer(),
                            process_mds_conn=self.mds_connection_initializer(),
                        )
                    )
                },
                tokamak=self.get_tokamak(),
                logger=self.logger,
            )
            shot_retriever.run(
                shot_creator_f=self._get_shot_data,
                shot_ids_list=shot_ids_list,
                await_complete=True,
            )
        else:
            shot_manager = self.get_shot_manager_cls()(
                process_database=self.database,
                process_mds_conn=self.mds_connection,
            )
            for shot_id in shot_ids_list:
                shot_data = self._get_shot_data(
                    shot_id=shot_id,
                    shot_manager=shot_manager,
                    shot_settings=shot_settings,
                )
                if shot_data is None:
                    self.logger.warning(
                        f"Not outputting data for shot {shot_id} due, data is None."
                    )
                else:
                    output_type_request.output_shot(
                        ResultOutputTypeRequestParams(
                            shot_id=shot_id,
                            result=shot_data,
                            database=self.database,
                            tokamak=self.get_tokamak(),
                            logger=self.logger,
                        )
                    )

        finish_output_type_request_params = FinishOutputTypeRequestParams(
            tokamak=self.get_tokamak(), logger=self.logger
        )
        results = output_type_request.get_results(finish_output_type_request_params)
        output_type_request.stream_output_cleanup(finish_output_type_request_params)
        return results
