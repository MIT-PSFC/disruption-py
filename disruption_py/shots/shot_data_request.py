from dataclasses import dataclass
from logging import Logger
import pandas as pd
import numpy as np
from disruption_py.shots.shot import Shot
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.method_caching import get_cached_method_params, is_cached_method, parameter_cached_method

from abc import ABC
from typing import Any, Callable, List


@dataclass
class ShotDataRequestParams:
    """Params passed by disruption_py to decorated methods.

    Attributes
    ----------
    shot : Any
		A reference to the shot object retrieving data.
    tokamak : Tokemak
        The tokemak for which the set times request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    shot : Shot
    tokamak : Tokamak
    logger : Logger
    
class ShotDataRequest(ABC):

    def setup(self, shot_data_request_params : ShotDataRequestParams):
        """Do any setup for the shot such as giving """
        pass
        
    def get_request_methods_for_tokamak(self, tokamak: Tokamak) -> List[Callable]:
        """Method used to determine which methods should be considered for execution given the tokamak.

        The default implementation returns all methods that have had the provided tokamak included in the tokamak parameter
        of there `cached_method` or `parameter_cached_method` decorator. This method may be overridden by subclasses if an
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
        return pd.DataFrame({"kappa_area": interp1(times, area/(np.pi * aminor**2), params.shot.times)})

