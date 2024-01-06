from dataclasses import dataclass
from logging import Logger
from disruption_py.shots.helpers.cached_method_params import is_cached_method
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.shots.helpers.cached_method_params import get_cached_method_params

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
        The tokamak for which the set times request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """
    shot_props : ShotProps
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
        for method_name in dir(self):
            method = getattr(self, method_name, None)
            if method is None or not is_cached_method(method):
                continue
            cached_method_params = get_cached_method_params(method, should_throw=True)
            if (cached_method_params.tokamaks is None or
                tokamak is cached_method_params.tokamaks or
                tokamak in cached_method_params.tokamaks):
                request_methods.append(method_name)
        return request_methods
