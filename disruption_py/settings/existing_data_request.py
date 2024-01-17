from abc import ABC, abstractmethod
from dataclasses import dataclass
import pandas as pd
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from typing import Dict, Union, Type
from logging import Logger

@dataclass
class ExistingDataRequestParams:
    """Params passed by disruption_py to _get_existing_data() method.

    Attributes
    ----------
    shot_id : int
        Shot Id for which to get existing data. Defaults to logbook.
    database : ShotDatabase
        Database object to use for getting existing data.
        A different database connection is used by each thread/process.
    tokamak : Tokemak
        The tokamak for which the data request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    shot_id : int
    database : ShotDatabase
    tokamak : Tokamak
    logger : Logger

ExistingDataRequestType = Union['ExistingDataRequest', str, pd.DataFrame, Dict[Tokamak, 'ExistingDataRequestType']]

class ExistingDataRequest(ABC):
    """ExistingDataRequest abstract class that should be inherited by all existing data request classes.
    """
    
    def get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_existing_data(params)
    
    @abstractmethod
    def _get_existing_data(self, params : ExistingDataRequestParams) -> pd.DataFrame:
        """Abstract method implemented by subclasses to get existing data for a given request as a pandas dataframe.
        
        Parameters
        ----------
        params : ExistingDataRequestParams
            Params that can be used to determine and retrieve existing data.
        
        Returns
        -------
        pd.DataFrame
            Pandas dataframe containing existing data.
        """
        pass
    
class ExistingDataRequestDict(ExistingDataRequest):
    def __init__(self, existing_data_request_dict : Dict[Tokamak, ExistingDataRequestType]):
        resolved_existing_data_request_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_existing_data_request(individual_request) 
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

# --8<-- [start:existing_data_request_dict]
_existing_data_request_mappings: Dict[str, ExistingDataRequest] = {
    "sql" : SQLExistingDataRequest(),
}
# --8<-- [end:existing_data_request_dict]

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
