
from typing import List, Callable
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
from disruption_py.handlers.multiprocessing_helper import MultiprocessingShotRetriever
from disruption_py.databases import CModDatabase
from disruption_py.shots import CModShot
import pandas as pd
import logging

class CModHandler:
    logger = logging.getLogger('disruption_py')
    
    def __init__(self, database_initializer : Callable[..., CModDatabase] = None, **kwargs):
        self.database_initializer = database_initializer
        if self.database_initializer is None:
            self.database_initializer = CModDatabase.default

    def get_shot_data(shot_id, sql_database=None, use_sql_table=True, **kwargs) -> pd.DataFrame:
        """
        Get shot data from CMOD.
        """
        class_logger = CModHandler.logger
        class_logger.info(f"starting {shot_id}")
        if use_sql_table:
            class_logger.info(f"retrieving sql data for {shot_id}")
            sql_shot_data = sql_database.get_shot_data(shot_ids=[shot_id])
        else:
            sql_shot_data = None
        disruption_time=sql_database.get_disruption_time(shot_id)
        try:
            shot = CModShot(shot_id=shot_id, existing_data=sql_shot_data, disruption_time=disruption_time, **kwargs)
            class_logger.info(f"completed {shot_id}")
            return shot.data
        except Exception as e:
            class_logger.error(f"failed {shot_id} with error {e}")
        return None

    def get_shots_data(self, shot_ids, multiprocessing=True, use_sql_table=True, **kwargs) -> List[pd.DataFrame]:
        """
        Get shot data from CMOD.
        """
        results = []
        if multiprocessing:
            shot_retriever = MultiprocessingShotRetriever(num_processes = 4, database_initializer_f=self.database_initializer)
            results = shot_retriever.run(
                shot_creator_f=CModHandler.get_shot_data, 
                shot_id_list=shot_ids,
                shot_args_dict={**kwargs, "use_sql_table": use_sql_table}, 
                should_finish=True
            )
        else:
            for shot_num in shot_ids:
                database = self.database_initializer()
                shot_data = self.get_shot_data(shot_num, sql_database=database, use_sql_table=use_sql_table, **kwargs)
                results.append(shot_data)
        return results
     