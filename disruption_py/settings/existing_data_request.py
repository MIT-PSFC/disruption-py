"""
Classes
----------
ExistingDataRequestParams
    Dataclass to hold parameters for a data request.
ExistingDataRequest
    Abstract class that must be subclassed by existing data request classes.
SQLExistingDataRequest
    Implementation for using the SQL database as the existing data request.
DFExistingDataRequest
    Implementation for passing existing data as a pandas DataFrame.

Attributes
----------
_existing_data_request_mappings (Dict[str, ExistingDataRequest]):
    Contains string mappings for the built-in existing data request classes
    Currently contains:
    "sql" : SQLExistingDataRequest

See Also:
    disruption_py.settings.shot_settings.ShotSettings
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
import pandas as pd
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokamak import Tokemak
from typing import Dict, Union, Type
from logging import Logger

@dataclass
class ExistingDataRequestParams:
    """Params passed by disruption_py to _get_existing_data() method.

    Attributes
    ----------
    shot_id : str
        Shot Id for which to get existing data.
    database : ShotDatabase
        Database object to use for getting existing data.
        A different database connection is used by each process.
    tokamak : Tokemak
        The tokemak for which the data request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    shot_id : str
    database : ShotDatabase
    tokamak : Tokemak
    logger : Logger

ExistingDataRequestType = Union['ExistingDataRequest', str, pd.DataFrame, Dict[Tokemak, 'ExistingDataRequestType']]

class ExistingDataRequest(ABC):
    """ExistingDataRequest abstract class that should be inherited by all existing data request classes.
    
    Methods
    -------
    _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame
        Abstract method implemented by subclasses to get existing data for a given request as a pandas dataframe.
    """
    
    def get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_existing_data(params)
    
    @abstractmethod
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        """Abstract method implemented by subclasses to get existing data for a given request as a pandas dataframe.
        
        Attributes
        ----------
        params : ExistingDataRequestParams
            Params that can be used to determine and retrieve existing data.
        """
        pass
    
class ExistingDataRequestDict(ExistingDataRequest):
    def __init__(self, existing_data_request_dict : Dict[Tokemak, ExistingDataRequestType]):
        resolved_existing_data_request_dict = {
            tokamak: resolve_existing_data_request(individual_request) 
            for tokamak, individual_request in existing_data_request_dict.items()
        }
        self.resolved_existing_data_request_dict = resolved_existing_data_request_dict
        
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        chosen_request = self.resolved_existing_data_request_dict.get(params.tokamak, None)
        if chosen_request is not None:
            return chosen_request.get_existing_data(params)
        else:
            params.logger.warning(f'No existing data request for tokamak {params.tokamak}')
            return None


class SQLExistingDataRequest(ExistingDataRequest):
    """Existing data request for retrieving data from SQL database."""
    
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        params.logger.info(f"retrieving sql data for {params.shot_id}")
        return params.database.get_shot_data(shot_ids=[params.shot_id])
    
class DFExistingDataRequest(ExistingDataRequest):
    """Existing data request for retrieving data from a pandas dataframe.
    
    Parameters
    ----------
    existing_data : pd.DataFrame
        The dataframe to use as the existing data.
    """

    def __init__(self, existing_data: pd.DataFrame):
        self.existing_data = existing_data
        
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        return self.existing_data
    
    
_existing_data_request_mappings: Dict[str, ExistingDataRequest] = {
    # do not include dataframe as requires an argument
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
