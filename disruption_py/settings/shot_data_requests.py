import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
from typing import List, Callable, Any
from logging import Logger
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.method_caching import parameter_cached_method

@dataclass
class ShotDataRequestParams:
    shot : Any
    existing_data : pd.DataFrame
    tokamak : Tokamak
    logger : Logger
    
class ShotDataRequest(ABC):
    
    def get_request_methods_for_tokamak(self, tokamak: Tokamak) -> List[Callable]:
        request_methods = []
        for method in dir(self):
            if hasattr(method, "tokamaks") and (tokamak in method.tokamaks or tokamak is method.tokamaks):
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
