from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokemak import Tokemak
from logging import Logger

@dataclass
class ShotNumberRequestParams:
    database : ShotDatabase
    tokemak : Tokemak
    logger : Logger
    
class ShotNumberRequest(ABC):
    
    def get_shot_numbers(self, params : ShotNumberRequestParams) -> List:
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._get_shot_numbers(params)
    
    @abstractmethod
    def does_use_database(self):
        pass
    
    @abstractmethod
    def _get_shot_numbers(self, params : ShotNumberRequestParams) -> List:
        """
        Default implementation of get_shot_numbers.
        Used for any tokemak not in tokemak_overrides.
        """
        pass

class ListShotNumberRequest(ShotNumberRequest):
    def __init__(self, shot_numbers):
        self.shot_numbers = shot_numbers

    def does_use_database(self):
        return False
    
    def _get_shot_numbers(self, params : ShotNumberRequestParams) -> List:
        return self.shot_numbers
    
class DatabaseShotNumberRequest(ShotNumberRequest):
    def __init__(self, sql_query, use_pandas=True):
        self.sql_query = sql_query
        self.use_pandas = use_pandas

    def does_use_database(self):
        return True
    
    def _get_shot_numbers(self, params : ShotNumberRequestParams) -> List:
        if self.use_pandas:
            query_result_df = params.database.query(query=self.sql_query, use_pandas=True)
            return query_result_df.iloc[:, 0].tolist()
        else:
            query_result = params.database.query(query=self.sql_query, use_pandas=False)
            return [row[0] for row in query_result]
     
_get_shot_numbers_request_mappings: Dict[str, ShotNumberRequest] = {
    # do not include classes that require initialization arguments
}

def shot_numbers_request_runner(shot_number_request, params : ShotNumberRequestParams):
    if isinstance(shot_number_request, ShotNumberRequest):
        return shot_number_request.get_shot_numbers(params)
    
    if isinstance(shot_number_request, int):
        return [shot_number_request]
    
    if isinstance(shot_number_request, str):
        timebase_request_object = _get_shot_numbers_request_mappings.get(shot_number_request, None)
        if timebase_request_object is not None:
            return timebase_request_object.get_times(params)
        
    if isinstance(shot_number_request, dict):
        chosen_request = shot_number_request.get(params.tokemak, None)
        if chosen_request is not None:
            return shot_numbers_request_runner(chosen_request, params)
        
    if isinstance(shot_number_request, list):
        all_results = []
        for request in shot_number_request:
            sub_result = shot_numbers_request_runner(request, params)
            if sub_result is not None:
                all_results.append(sub_result)
        
        return [shot_num for sub_list in all_results for shot_num in sub_list]
    
    raise ValueError("Invalid shot number request")
