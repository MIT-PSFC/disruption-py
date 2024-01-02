from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Union, Type
import pandas as pd
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from logging import Logger

import disruption_py.data
from disruption_py.utils.utils import without_duplicates
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources


@dataclass
class ShotIdsRequestParams:
    """Params passed by disruption_py to _get_shot_ids() method.

    Attributes
    ----------
    database : ShotDatabase
        Database object to use for getting shot ids. 
        A different database connection is used by each process.
        Defaults to logbook.
    tokamak : Tokemak
        The tokemak for which the data request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    database : ShotDatabase
    tokamak : Tokamak
    logger : Logger

class ShotIdsRequest(ABC):
    """ShotIdsRequest abstract class that should be inherited by all shot id request classes."""
    
    def get_shot_ids(self, params : ShotIdsRequestParams) -> List:
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_shot_ids(params)
    
    @abstractmethod
    def does_use_database(self):
        """Abstract method implemented by subclasses to signify whether a database connection is used.
        This allows disruption_py to only initialize a database connection in the parent process if 
        it is being used.
        """
    
    @abstractmethod
    def _get_shot_ids(self, params : ShotIdsRequestParams) -> List:
        """Abstract method implemented by subclasses to get shot ids for the given request params as a list.
        
        Parameters
        ----------
        params : ShotIdsRequestParams
            Params that can be used to determine shot ids.
        """
        pass

class IncludedShotIdsRequest(ShotIdsRequest):
    """Use the shot IDs from one of the provided data files.
    
    Directly passing a key from the _get_shot_ids_request_mappings dictionary as a string will
    automatically crete a new IncludedShotIdsRequest object with that data_file_name.
    
    Parameters
    ----------
    data_file_name : str
        The name of the datafile that should be used to retrieve shot_ids.
    """
    
    def __init__(self, data_file_name):
        with importlib_resources.path(disruption_py.data, data_file_name) as p:
            data_file_name = str(p)
        self.shot_ids = pd.read_csv(data_file_name, header=None).iloc[:, 0].values.tolist()

    def does_use_database(self):
        return False
    
    def _get_shot_ids(self, params : ShotIdsRequestParams) -> List:
        return self.shot_ids 

class FileShotIdsRequest(ShotIdsRequest):
    """Use a list of shot IDs from the provided file path, this may be any file readable by pandas read_csv.
    
    Directly passing a file path as a string to the shot id request with the file name suffixed by txt or csv
    will automatically create a new FileShotIdsRequest object with that file path.
    
    Parameters
    ----------
    file_path : str
        The file path of the file that should be used for retrieving shot ids.
    column_index : int
        The index of the column that should be read. For text files, this should be 0. Defaults to 0.
    """
    
    def __init__(self, file_path, column_index=0):
        self.shot_ids = pd.read_csv(file_path, header=None).iloc[:, column_index].values.tolist()

    def does_use_database(self):
        return False
    
    def _get_shot_ids(self, params : ShotIdsRequestParams) -> List:
        return self.shot_ids 
    
class DatabaseShotIdsRequest(ShotIdsRequest):
    """Use an sql query of the database to retrieve the shot ids.
    
    Parameters
    ----------
    sql_query : str
        The sql query that should be used for retrieving shot ids.
    use_pandas : bool
        Whether Pandas should be used to do the sql query. Defaults to true.
    """

    def __init__(self, sql_query, use_pandas=True):
        self.sql_query = sql_query
        self.use_pandas = use_pandas

    def does_use_database(self):
        return True
    
    def _get_shot_ids(self, params : ShotIdsRequestParams) -> List:
        if self.use_pandas:
            query_result_df = params.database.query(query=self.sql_query, use_pandas=True)
            return query_result_df.iloc[:, 0].tolist()
        else:
            query_result = params.database.query(query=self.sql_query, use_pandas=False)
            return [row[0] for row in query_result]

# --8<-- [start:get_shot_ids_request_dict]
_get_shot_ids_request_mappings: Dict[str, ShotIdsRequest] = {
    "paper": IncludedShotIdsRequest("paper_shotlist.txt"),
    "disr": IncludedShotIdsRequest("train_disr.txt"),
    "nondisr": IncludedShotIdsRequest("train_nondisr.txt"),
}
# --8<-- [end:get_shot_ids_request_dict]

# --8<-- [start:file_suffix_to_shot_ids_request_dict]
_file_suffix_to_shot_ids_request : Dict[str, Type[ShotIdsRequest]] = {
    ".txt" : FileShotIdsRequest,
    ".csv" : FileShotIdsRequest,
}
# --8<-- [end:file_suffix_to_shot_ids_request_dict]

ShotIdsRequestType = Union['ShotIdsRequest', int, str, Dict[Tokamak, 'ShotIdsRequestType'], List['ShotIdsRequestType']]

def shot_ids_request_runner(shot_ids_request, params : ShotIdsRequestParams):
    if isinstance(shot_ids_request, ShotIdsRequest):
        shot_ids = shot_ids_request.get_shot_ids(params)
    
    elif isinstance(shot_ids_request, int) or (isinstance(shot_ids_request, str) and shot_ids_request.isdigit()):
        shot_ids = [shot_ids_request]
    
    elif isinstance(shot_ids_request, str):
        shot_ids_request_object = _get_shot_ids_request_mappings.get(shot_ids_request, None)
        if shot_ids_request_object is not None:
            shot_ids = shot_ids_request_object.get_shot_ids(params)
        
    elif isinstance(shot_ids_request, str):
        # assume that it is a file path
       for suffix, shot_ids_request_type in _file_suffix_to_shot_ids_request.items():
           if shot_ids_request.endswith(suffix):
               shot_ids = shot_ids_request_type(shot_ids_request).get_shot_ids(params)
        
    elif isinstance(shot_ids_request, dict):
        shot_ids_request = {
            map_string_to_enum(tokamak, Tokamak): shot_ids_request_mapping 
            for tokamak, shot_ids_request_mapping in shot_ids_request.items()
        }
        chosen_request = shot_ids_request.get(params.tokamak, None)
        if chosen_request is not None:
            shot_ids = shot_ids_request_runner(chosen_request, params)
        
    elif isinstance(shot_ids_request, list):
        all_results = []
        for request in shot_ids_request:
            sub_result = shot_ids_request_runner(request, params)
            if sub_result is not None:
                all_results.append(sub_result)
        
        shot_ids = [shot_id for sub_list in all_results for shot_id in sub_list]
    else:
        raise ValueError("Invalid shot id request")
    
    return without_duplicates(shot_ids)
