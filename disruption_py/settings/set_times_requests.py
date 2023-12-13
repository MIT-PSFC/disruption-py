from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, Callable, Union
from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.utils.mappings.tokamak import Tokemak
from logging import Logger

@dataclass
class SetTimesRequestParams:
    tree_manager : TreeManager
    tokamak : Tokemak
    logger : Logger
    disruption_time : float = None

SetTimesRequestType = Union['SetTimesRequest', str, Dict[Tokemak, 'SetTimesRequestType']]

class SetTimesRequest(ABC):
    '''
    Represents a request for setting the times of a timebase for a shot
    
    Set tokamak_overrides to override the default behaviour for a tokamak.
    Subclasses must implement _get_times for the default case.

    note: Object should not be modified after being passed to a handler, 
    to avoid issues with shared memory when multiprocessing.
    '''
            
    def get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        if hasattr(self, 'tokamak_overrides'):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_times(params)
    
    @abstractmethod
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        """
        Default implementation of get_timebase.
        Used for any tokamak not in tokamak_overrides.
        """
        pass
    
class SetTimesRequestDict(SetTimesRequest):
    def __init__(self, set_times_request_dict : Dict[Tokemak, SetTimesRequestType]):
        resolved_set_times_request_dict = {
            tokamak: resolve_set_times_request(individual_request) 
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
    def __init__(self, times):
        self.times = times

    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        return self.times
    
class EfitSetTimesRequest(SetTimesRequest):
    def __init__(self):
        self.tokamak_overrides = {
            Tokemak.CMOD: self.cmod_times
        }

    def cmod_times(self, params : SetTimesRequestParams):
        efit_tree = params.tree_manager.tree_from_nickname("efit_tree")
        efit_tree_name = params.tree_manager.tree_name_of_nickname("efit_tree")
        if efit_tree_name == 'analysis':
            try:
                return efit_tree.getNode(r"\analysis::efit_aeqdsk:time").getData().data().astype('float64', copy=False)
            except Exception as e:
                return efit_tree.getNode(r"\analysis::efit:results:a_eqdsk:time").getData().data().astype('float64', copy=False)
        else:
            return efit_tree.getNode(fr"\{efit_tree_name}::top.results.a_eqdsk:time").getData().data().astype('float64', copy=False)
            
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        raise ValueError("EFIT timebase request not implemented")
    
class MagneticsSetTimesRequest(SetTimesRequest):
    
    def __init__(self, timestep = 0.004):
        self.timestep = timestep
        
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        try:
            magnetics_tree = params.tree_manager.open_tree(tree_name='magnetics')
            magnetic_times = magnetics_tree.getNode('\ip').dim_of().data().astype('float64', copy=False)
            # interpolate self._times to be every [timestep] seconds
            return np.arange(magnetic_times[0], magnetic_times[-1], self.timestep).astype('float64', copy=False)
        except:
            params.logger.info("Failed to set up magnetic timebase")
            return None
        

_set_times_request_mappings: Dict[str, SetTimesRequest] = {
    # do not include times list as requires an argument
    "efit" : EfitSetTimesRequest(),
    "magnetics004" : MagneticsSetTimesRequest(timestep = 0.004),
}

def resolve_set_times_request(set_times_request : SetTimesRequestType) -> SetTimesRequest:
    if isinstance(set_times_request, SetTimesRequest):
        return set_times_request
    
    if isinstance(set_times_request, str):
        set_times_request_object = _set_times_request_mappings.get(set_times_request, None)
        if set_times_request_object is not None:
            return set_times_request_object
        
    if isinstance(set_times_request, np.ndarray):
        return ListSetTimesRequest(set_times_request)
    
    if isinstance(set_times_request, pd.Series):
        return ListSetTimesRequest(set_times_request.to_numpy())
        
    if isinstance(set_times_request, dict):
        return SetTimesRequestDict(set_times_request)
    
    raise ValueError("Invalid set times request")