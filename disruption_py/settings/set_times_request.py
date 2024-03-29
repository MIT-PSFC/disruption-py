from abc import ABC, abstractmethod
import traceback
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, Union
from disruption_py.databases.database import ShotDatabase
from disruption_py.mdsplus_integration.mds_connection import MDSConnection
from disruption_py.utils.constants import MAX_SHOT_TIME
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak
from logging import Logger

@dataclass
class SetTimesRequestParams:
    """Params passed by disruption_py to _get_times() method.

    Attributes
    ----------
    shot_id : int
        The shot id for the timebase being created
    mds_conn : MDSConnection
        The connection to MDSplus, that can be used to retrieve MDSplus data.
    existing_data : pd.DataFrame
        Pre-filled data given to disruption_py.
    database : ShotDatabase
        Database object with connection to sql database.
    disruption_time : float
        The time when the shot disrupted or None if no disruption occured.
    tokamak : Tokemak
        The tokamak for which the set times request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    shot_id : int
    mds_conn : MDSConnection
    existing_data : pd.DataFrame
    database : ShotDatabase
    disruption_time : float
    tokamak : Tokamak
    logger : Logger

SetTimesRequestType = Union['SetTimesRequest', str, Dict[Tokamak, 'SetTimesRequestType']]

class SetTimesRequest(ABC):
    """SetTimesRequest abstract class that should be inherited by all set times request classes.
    """
            
    def get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_times(params)
    
    @abstractmethod
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        """Abstract method implemented by subclasses to get timebase as list.
        The timebase can be set to be automatically restricted to a subdomain of the
        provided times via the signal_domain argument in the ShotSettings object.
        
        Parameters
        ----------
        params : SetTimesRequestParams
            Params that can be used to determine and retrieve the timebase.
        
        Returns
        -------
        np.ndarray
            Numpy array containing times in the timebase.
        """
        pass
    
class SetTimesRequestDict(SetTimesRequest):
    """
    Utility class that is automatically used when a dicationary is passed as the `set_times_request` parameter in `ShotSettings`.
    
    Parameters
    ----------
    set_times_request_dict : dict[Tokamak, SetTimesRequestType]
        A dictionary mapping tokamak type strings to the desired set times request for that tokamak.  E.g. `{'cmod': 'efit'}`.
        Any other option passable to the `set_times_request` parameter in `ShotSettings` may be used.
    """
    def __init__(self, set_times_request_dict : Dict[Tokamak, SetTimesRequestType]):
        resolved_set_times_request_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_set_times_request(individual_request) 
            for tokamak, individual_request in set_times_request_dict.items()
        }
        self.resolved_set_times_request_dict = resolved_set_times_request_dict
        
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        chosen_request = self.resolved_set_times_request_dict.get(params.tokamak, None)
        if chosen_request is not None:
            return chosen_request.get_times(params)
        else:
            params.logger.warning(f'No get times request for tokamak {params.tokamak}')
            return None

class ListSetTimesRequest(SetTimesRequest):
    """ 
    A list of times to use as the timebase. 
    
    A utility class that is automatically used when a list, numpy array, or pandas series is
    passed as the `set_times_request` parameter in `ShotSettings`.
    """
     
    def __init__(self, times):
        self.times = times

    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        return self.times

class ExistingDataSetTimesRequest(SetTimesRequest):
    """
    Get times request for using the existing data timebase.
    
    If no existing data exists for the shot, then the fallback_set_times_request is used.
    """
    def __init__(self, fallback_set_times_request : SetTimesRequestType) -> None:
        self.fallback_set_times_request = resolve_set_times_request(fallback_set_times_request)

    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        if params.existing_data is not None:
            # set timebase to be the timebase of existing data
            try:
                times = params.existing_data['time'].to_numpy()
                # Check if the timebase is in ms instead of s
                if times[-1] > MAX_SHOT_TIME:
                    times /= 1000  # [ms] -> [s]
                return times
            except KeyError as e:
                params.logger.warning(
                    f"[Shot {params.shot_id}]: Shot constructor was passed data but no timebase.")
                params.logger.debug(
                    f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        else:
            return self.fallback_set_times_request.get_times(params)
        
        
class EfitSetTimesRequest(SetTimesRequest):
    """ Get times request for using the EFIT timebase.  """
    
    def __init__(self):
        self.tokamak_overrides = {
            Tokamak.CMOD: self.cmod_times
        }

    def cmod_times(self, params : SetTimesRequestParams):
        efit_tree_name = params.mds_conn.get_tree_name_of_nickname("_efit_tree")
        if efit_tree_name == 'analysis':
            try:
                return params.mds_conn.get_data(r"\analysis::efit_aeqdsk:time", tree_name="_efit_tree", astype="float64")
            except Exception as e:
                return params.mds_conn.get_data(r"\analysis::efit:results:a_eqdsk:time", tree_name="_efit_tree", astype="float64")
        else:
            return params.mds_conn.get_data(fr"\{efit_tree_name}::top.results.a_eqdsk:time", tree_name="_efit_tree", astype="float64")
            
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        raise ValueError("EFIT timebase request not implemented")
    
        
class SignalSetTimesRequest(SetTimesRequest):
    def __init__(self, tree_name : str, signal_path : str):
        self.tree_name = tree_name
        self.signal_path = signal_path
        
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        try:
            signal_time, = params.mds_conn.get_dims(self.signal_path, tree_name=self.tree_name)
            return signal_time.astype('float64', copy=False)
        except Exception as e:
            params.logger.error(f"Failed to set up timebase for signal {self.signal_path}")
            raise Exception(e)
        

# --8<-- [start:set_times_request_dict]
_set_times_request_mappings: Dict[str, SetTimesRequest] = {
    "efit" : EfitSetTimesRequest(),
}
# --8<-- [end:set_times_request_dict]

def resolve_set_times_request(set_times_request : SetTimesRequestType) -> SetTimesRequest:
    if isinstance(set_times_request, SetTimesRequest):
        return set_times_request
    
    if isinstance(set_times_request, str):
        set_times_request_object = _set_times_request_mappings.get(set_times_request, None)
        if set_times_request_object is not None:
            return set_times_request_object
        
    if isinstance(set_times_request, np.ndarray) or isinstance(set_times_request, list):
        return ListSetTimesRequest(set_times_request)
    
    if isinstance(set_times_request, pd.Series):
        return ListSetTimesRequest(set_times_request.to_numpy())
        
    if isinstance(set_times_request, dict):
        return SetTimesRequestDict(set_times_request)
    
    raise ValueError("Invalid set times request")