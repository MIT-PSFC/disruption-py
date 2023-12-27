from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, Union
from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak
from logging import Logger

@dataclass
class SetTimesRequestParams:
    """Params passed by disruption_py to _get_times() method.

    Attributes
    ----------
    tree_manager : TreeManager
        Tree manager which can be used to retrieve data from MDSplus.
    database : ShotDatabase
        Database object to use for getting timebase from sql database.
        A different database connection is used by each process.
    tokamak : Tokemak
        The tokemak for which the set times request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    tree_manager : TreeManager
    tokamak : Tokamak
    logger : Logger
    disruption_time : float = None

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
        
        Attributes
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
    def __init__(self, times):
        self.times = times

    def _get_times(self, params : SetTimesRequestParams) -> np.ndarray:
        return self.times
    
class EfitSetTimesRequest(SetTimesRequest):
    """Get times request for using the EFIT timebase.
    """
    
    def __init__(self):
        self.tokamak_overrides = {
            Tokamak.CMOD: self.cmod_times
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
    """Get times request for using the start and end times of the magnetics tree, with a 
    custom timestep passed in the constructor.
    
    Parameters
    ----------
    timestep : float
        The timestep to use for the magnetics timebase.
    """
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
        

# --8<-- [start:set_times_request_dict]
_set_times_request_mappings: Dict[str, SetTimesRequest] = {
    "efit" : EfitSetTimesRequest(),
    "magnetics004" : MagneticsSetTimesRequest(timestep = 0.004),
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