
from typing import Callable, Any
import traceback
from disruption_py.handlers.multiprocessing_helper import MultiprocessingShotRetriever
from disruption_py.settings.shot_id_requests import ShotIdRequestParams, ShotIdRequestType, shot_ids_request_runner
from disruption_py.settings.existing_data_request import ExistingDataRequest, ExistingDataRequestParams
from disruption_py.settings.output_type_requests import ResultOutputTypeRequestParams, FinishOutputTypeRequestParams
from disruption_py.settings import ShotDataRequestParams, ShotSettings
from disruption_py.utils.mappings.tokamak import Tokemak
from disruption_py.databases import CModDatabase
from disruption_py.shots import CModShot
from disruption_py.shots.populate_shot import populate_shot
import pandas as pd
import logging

class CModHandler:
    """
    Brief description of the class.

    Parameters
    ----------
    database_initializer : Callable[..., CModDatabase]
        When run returns a new database object for the handler. The function must create a new database 
        connection instead of reusing an existing one, as the handler may initalize multiple connections 
        across different processes. Defaults to CModDatabase.default.

    Attributes
    ----------
    logger : Logger
        The logger used for disruption_py.

    Methods
    -------
    get_shot_data(shot_id, sql_database=None, shot_settings)
        Static method used to get data for a single shot from CMOD. May be run across different processes.
    get_shots_data(shot_id_request, shot_settings, num_processes)
        Instance method used to get shot data for all shots from shot_id_request from CMOD.
    """
    logger = logging.getLogger('disruption_py')
    
    def __init__(self, database_initializer : Callable[..., CModDatabase] = None, **kwargs):
        self.database_initializer = database_initializer
        if self.database_initializer is None:
            self.database_initializer = CModDatabase.default

    @property
    def database(self):
        """
        Accessor for the database.
        """
        if not hasattr(self, '_database'):
            self._database = self.database_initializer()
        return self._database
    
    @staticmethod
    def get_shot_data(shot_id, sql_database=None, shot_settings: ShotSettings=None) -> pd.DataFrame:
        """
        Get data for a single shot from CMOD. May be run across different processes.
        """
        tokamak = Tokemak.CMOD
        class_logger = CModHandler.logger
        class_logger.info(f"starting {shot_id}")
        if shot_settings.existing_data_request is not None:
            existing_data_request_params = ExistingDataRequestParams(
                shot_id=str(shot_id),
                database=sql_database,
                tokamak=tokamak, 
                logger=class_logger,
            )
            existing_data = shot_settings.existing_data_request.get_existing_data(existing_data_request_params)
        else:
            existing_data = None
        disruption_time=sql_database.get_disruption_time(shot_id)
        try:
            shot = CModShot(shot_id=shot_id, existing_data=existing_data, disruption_time=disruption_time, shot_settings=shot_settings)
            retrieved_data = populate_shot(shot_run_settings=shot_settings, params=ShotDataRequestParams(shot, existing_data, tokamak, class_logger))
            shot.cleanup()
            class_logger.info(f"completed {shot_id}")
            return retrieved_data
        except Exception as e:
            class_logger.warning(f"[Shot {shot_id}]: fatal error {traceback.format_exc()}")
            class_logger.error(f"failed {shot_id} with error {e}")
        return None

    def get_shots_data(
        self,
        shot_id_request : ShotIdRequestType,
        shot_settings : ShotSettings = None,
        num_processes: int = 1,
    ) -> Any:
        """
        Get shot data for all shots from shot_id_request from CMOD.
        
        Attributes
        ----------
        shot_id_request : ShotIdRequestType
            Data retrieved for all shot_ids specified by the request. See ShotIdRequest for more details.
        shot_settings : ShotSettings
            The settings that each shot uses when retrieving data. See ShotSettings for more details.
            If None, the default values of each setting in ShotSettings is used.
        num_processes : int
            The number of processes to use for data retrieval. If 1, the data is retrieved in serial. 
            If > 1, the data is retrieved in parallel.
            
        Returns
        -------
        Any
            The value of OutputTypeRequest.get_results, where OutputTypeRequest is specified in 
            shot_settings. See OutputTypeRequest for more details.
        """
        tokamak = Tokemak.CMOD
        
        if shot_settings is None:
            shot_settings = ShotSettings()
        shot_settings.resolve()
        
        shot_id_request_params = ShotIdRequestParams(self.database, tokamak, self.logger)
        shot_id_list = shot_ids_request_runner(shot_id_request, shot_id_request_params)
        
        if num_processes > 1:
            shot_retriever = MultiprocessingShotRetriever(
                database_initializer_f=self.database_initializer,
                num_processes=num_processes,
                shot_settings=shot_settings,
                tokamak = tokamak,
                logger = self.logger,
            )
            results = shot_retriever.run(
                shot_creator_f=CModHandler.get_shot_data, 
                shot_id_list=shot_id_list,
            )
        else:
            for shot_id in shot_id_list:
                shot_data = CModHandler.get_shot_data(
                    shot_id=shot_id, 
                    sql_database=self.database, 
                    shot_settings=shot_settings
                )
                shot_settings.output_type_request.output_shot(ResultOutputTypeRequestParams(shot_data, Tokemak.CMOD, self.logger))
            
            finish_output_type_request_params = FinishOutputTypeRequestParams(tokamak, self.logger)
            shot_settings.output_type_request.stream_output_cleanup(finish_output_type_request_params)
            results = shot_settings.output_type_request.get_results(finish_output_type_request_params)
        return results
     