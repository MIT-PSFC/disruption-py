from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Union, Type
import pandas as pd
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokamak import Tokemak
from logging import Logger

import disruption_py.data
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources


@dataclass
class ShotIdRequestParams:
    database : ShotDatabase
    tokamak : Tokemak
    logger : Logger

class ShotIdRequest(ABC):
    
    def get_shot_ids(self, params : ShotIdRequestParams) -> List:
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_shot_ids(params)
    
    @abstractmethod
    def does_use_database(self):
        pass
    
    @abstractmethod
    def _get_shot_ids(self, params : ShotIdRequestParams) -> List:
        """
        Default implementation of get_shot_ids.
        Used for any tokamak not in tokamak_overrides.
        """
        pass

class IncludedShotIdRequest(ShotIdRequest):
    '''
    Use a list of shot IDs from one of the provided data files
    '''
    def __init__(self, data_file_name):
        with importlib_resources.path(disruption_py.data, data_file_name) as p:
            data_file_name = str(p)
        self.shot_ids = pd.read_csv(data_file_name, header=None).iloc[:, 0].values.tolist()

    def does_use_database(self):
        return False
    
    def _get_shot_ids(self, params : ShotIdRequestParams) -> List:
        return self.shot_ids 
    
class ListShotIdRequest(ShotIdRequest):
    def __init__(self, shot_ids):
        self.shot_ids = shot_ids

    def does_use_database(self):
        return False
    
    def _get_shot_ids(self, params : ShotIdRequestParams) -> List:
        return self.shot_ids

class FileShotIdRequest(ShotIdRequest):
    '''
    Use a list of shot IDs from the provided file name, this may be any file readable by pandas read_csv.
    '''
    def __init__(self, file_path, column_index=0):
        self.shot_ids = pd.read_csv(file_path, header=None).iloc[:, column_index].values.tolist()

    def does_use_database(self):
        return False
    
    def _get_shot_ids(self, params : ShotIdRequestParams) -> List:
        return self.shot_ids 
    
class DatabaseShotIdRequest(ShotIdRequest):
    def __init__(self, sql_query, use_pandas=True):
        self.sql_query = sql_query
        self.use_pandas = use_pandas

    def does_use_database(self):
        return True
    
    def _get_shot_ids(self, params : ShotIdRequestParams) -> List:
        if self.use_pandas:
            query_result_df = params.database.query(query=self.sql_query, use_pandas=True)
            return query_result_df.iloc[:, 0].tolist()
        else:
            query_result = params.database.query(query=self.sql_query, use_pandas=False)
            return [row[0] for row in query_result]
     
_get_shot_id_request_mappings: Dict[str, ShotIdRequest] = {
    # do not include classes that require initialization arguments
    "paper": IncludedShotIdRequest("paper_shotlist.txt"),
    "disr": IncludedShotIdRequest("train_disr.txt"),
    "nondisr": IncludedShotIdRequest("train_nondisr.txt"),
}

_file_suffix_to_shot_id_request : Dict[str, Type[ShotIdRequest]] = {
    ".txt" : FileShotIdRequest,
    ".csv" : FileShotIdRequest,
} 

ShotIdRequestType = Union['ShotIdRequest', int, str, Dict[Tokemak, 'ShotIdRequestType'], List['ShotIdRequestType']]

def shot_ids_request_runner(shot_id_request, params : ShotIdRequestParams):
    if isinstance(shot_id_request, ShotIdRequest):
        return shot_id_request.get_shot_ids(params)
    
    if isinstance(shot_id_request, int) or (isinstance(shot_id_request, str) and shot_id_request.isdigit()):
        return [shot_id_request]
    
    if isinstance(shot_id_request, str):
        shot_id_request_object = _get_shot_id_request_mappings.get(shot_id_request, None)
        if shot_id_request_object is not None:
            return shot_id_request_object.get_shot_ids(params)
        
    if isinstance(shot_id_request, str):
        # assume that it is a file path
       for suffix, shot_id_request_type in _file_suffix_to_shot_id_request.items():
           if shot_id_request.endswith(suffix):
               return shot_id_request_type(shot_id_request).get_shot_ids(params)
        
    if isinstance(shot_id_request, dict):
        chosen_request = shot_id_request.get(params.tokamak, None)
        if chosen_request is not None:
            return shot_ids_request_runner(chosen_request, params)
        
    if isinstance(shot_id_request, list):
        all_results = []
        for request in shot_id_request:
            sub_result = shot_ids_request_runner(request, params)
            if sub_result is not None:
                all_results.append(sub_result)
        
        return [shot_id for sub_list in all_results for shot_id in sub_list]
    
    raise ValueError("Invalid shot id request")
