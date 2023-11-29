from abc import ABC, abstractmethod
import numpy as np
from dataclasses import dataclass
from typing import Dict, Callable
from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.utils.mappings.tokemak import Tokemak
from logging import Logger

@dataclass
class SetTimesRequestParams:
    tree_manager : TreeManager
    tokemak : Tokemak
    logger : Logger
    disruption_time : float = None

class SetTimesRequest(ABC):
    '''
    Represents a request for setting the times of a taimebase for a shot
    
    Set tokemak_overrides to override the default behavior for a tokemak.
    Subclasses must implement _get_times for the default case.

    note: Object should not be modified after being passed to a handler, 
    to avoid issues with shared memory when multiprocessing.
    '''
            
    def get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._get_times(params)
    
    @abstractmethod
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        """
        Default implementation of get_timebase.
        Used for any tokemak not in tokemak_overrides.
        """
        pass

class ListSetTimesRequest(SetTimesRequest):
    def __init__(self, times):
        self.times = times

    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        return self.times
    
class EfitSetTimesRequest(SetTimesRequest):
    def __init__(self):
        self.tokemak_overrides = {
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
            return efit_tree.getNode(fr"\{efit_tree_name}::efit.results.a_eqdsk:time").getData().data().astype('float64', copy=False)
            
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        raise ValueError("EFIT timebase request not implemented")
    
class MagneticsSetTimesRequest(SetTimesRequest):
    
    def __init__(self, timestep = 0.004):
        self.timestep = timestep
        
    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        try:
            magnetics_tree = params.tree_manager.open_tree(tree_name='magnetics')
            magnetic_times = magnetics_tree.getNode('\ip').dim_of().data().astype('float64', copy=False)
            # # interpolate self._times to be every [timestep] seconds
            return np.arange(magnetic_times[0], magnetic_times[-1], self.timestep).astype('float64', copy=False)
        except:
            params.logger.info("Failed to set up magnetic timebase")
            return None
        

_set_times_request_mappings: Dict[str, SetTimesRequest] = {
    # do not include times list as requires an argument
    "efit" : EfitSetTimesRequest(),
    "magnetics" : MagneticsSetTimesRequest(),
}

def set_times_request_runner(set_times_request, params : SetTimesRequestParams) -> np.ndarray:
    if isinstance(set_times_request, SetTimesRequest):
        return set_times_request.get_times(params)
    
    if isinstance(set_times_request, str):
        set_times_request_object = _set_times_request_mappings.get(set_times_request, None)
        if set_times_request_object is not None:
            return set_times_request_object.get_times(params)
        
    if isinstance(set_times_request, dict):
        chosen_request = set_times_request.get(params.tokemak, None)
        if chosen_request is not None:
            return set_times_request_runner(chosen_request, params)
    
    raise ValueError("Invalid timebase request")