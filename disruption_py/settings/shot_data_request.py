import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
from typing import List, Callable, Any
from logging import Logger
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.method_caching import parameter_cached_method, is_cached_method, get_cached_method_params

@dataclass
class ShotDataRequestParams:
    """Params passed by disruption_py to decorated methods.

    Attributes
    ----------
    shot : Shot
        An instance of the tokamak's shot class for the shot that data is being retrieved for.
        This shot object should be used to access the shot's MDSplus trees and timebase using
        the `get_tree_manager()` and `get_times()` methods respectively.
    existing_data : pd.DataFrame
        Data provided to disruption_py for the given shot in the `existing_data_request` parameter of `shot_settings`.
    tokamak : Tokemak
        The tokemak for which the set times request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    shot : Any
    existing_data : pd.DataFrame
    tokamak : Tokamak
    logger : Logger
    
class ShotDataRequest(ABC):
    
    def get_request_methods_for_tokamak(self, tokamak: Tokamak) -> List[Callable]:
        """Method used to determine which methods should be considered for execution given the tokamak.
        
        The default implementation returns all methods that have had the provided tokamak included in the tokamak parameter
        of there `cached_method` or `parameter_cached_method` decorator. This method may ve overridden by subclasses if an
        alternate scheme of determining which methods can be executed for each tokamak is required.

        Parameters
        ----------
        tokamak : Tokamak
            The tokamak for which we are retrieving elligible methods.

        Returns
        -------
        List[Callable]
            A list of methods that can be considered for execution for the given tokamak.
        """
        request_methods = []
        for method in dir(self):
            if not is_cached_method:
                continue
            cached_method_params = get_cached_method_params(method, should_throw=True)
            if (cached_method_params.tokamaks is None or
                tokamak in cached_method_params.tokamaks or 
                tokamak is cached_method_params.tokamaks):
                request_methods.append(method)
        return request_methods
    
class KappaArea(ShotDataRequest):
    
    @parameter_cached_method(columns=["kappa_area"], used_trees=["efit_tree"], tokamaks=Tokamak.CMOD)
    def _get_kappa_area(self, params):
        aminor = params.shot.efit_tree.getNode(
            r'\efit_aeqdsk:aminor').getData().data().astype('float64', copy=False)
        area = params.shot.efit_tree.getNode(
            r'\efit_aeqdsk:area').getData().data().astype('float64', copy=False)
        times = params.shot.efit_tree.getNode(
            r'\efit_aeqdsk:time').getData().data().astype('float64', copy=False)

        aminor[aminor <= 0] = 0.001  # make sure aminor is not 0 or less than 0
        # make sure area is not 0 or less than 0
        area[area <= 0] = 3.14*0.001**2
        return pd.DataFrame({"kappa_area": interp1(times, area/(np.pi * aminor**2), params.shot.get_times())})
