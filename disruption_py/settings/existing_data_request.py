from abc import ABC, abstractmethod
from dataclasses import dataclass
import pandas as pd
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokemak import Tokemak
from typing import Dict, Union, Type
from logging import Logger

@dataclass
class ExistingDataRequestParams:
    shot_id : str
    database : ShotDatabase
    tokemak : Tokemak
    logger : Logger

ExistingDataRequestType = Union['ExistingDataRequest', str, pd.DataFrame, Dict[Tokemak, 'ExistingDataRequestType']]

class ExistingDataRequest(ABC):
    
    def get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._get_existing_data(params)
    
    @abstractmethod
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        pass
    
class ExistingDataRequestDict(ExistingDataRequest):
    def __init__(self, existing_data_request_dict : Dict[Tokemak, ExistingDataRequestType]):
        resolved_existing_data_request_dict = {
            tokemak: resolve_existing_data_request(individual_request) 
            for tokemak, individual_request in existing_data_request_dict.items()
        }
        self.resolved_existing_data_request_dict = resolved_existing_data_request_dict
        
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        chosen_request = self.resolved_existing_data_request_dict.get(params.tokemak, None)
        if chosen_request is not None:
            return chosen_request.get_existing_data(params)
        else:
            params.logger.warning(f'No existing data request for tokemak {params.tokemak}')
            return None


class SQLExistingDataRequest(ExistingDataRequest):
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        params.logger.info(f"retrieving sql data for {params.shot_id}")
        return params.database.get_shot_data(shot_ids=[params.shot_id])
    
class DFExistingDataRequest(ExistingDataRequest):
    def __init__(self, existing_data: pd.DataFrame):
        self.existing_data = existing_data
        
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        return self.existing_data
    
    
_existing_data_request_mappings: Dict[str, ExistingDataRequest] = {
    # do not include times list as requires an argument
    "sql" : SQLExistingDataRequest(),
}

def resolve_existing_data_request(existing_data_request : ExistingDataRequestType) -> ExistingDataRequest:
    if existing_data_request is None:
        return None
    
    if isinstance(existing_data_request, ExistingDataRequest):
        return existing_data_request
    
    if isinstance(existing_data_request, str):
        existing_data_request_object = _existing_data_request_mappings.get(existing_data_request, None)
        if existing_data_request_object is not None:
            return existing_data_request_object
        
    if isinstance(existing_data_request, pd.DataFrame):
        return DFExistingDataRequest(existing_data_request)
        
    if isinstance(existing_data_request, dict):
        return ExistingDataRequestDict(existing_data_request)
    
    raise ValueError("Invalid set times request")
