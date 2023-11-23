from dataclasses import dataclass
from abc import ABC, abstractmethod
from enum import Enum
import numpy as np

from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.utils.mappings.tokemak import Tokemak
from logging import Logger

@dataclass
class SignalDomainSubrequestParams:
    tree_manager : TreeManager
    tokemak : Tokemak
    logger : Logger
    kwargs : dict


class SignalDomainSubrequest(ABC):

    def get_signal_domain(self, params : SignalDomainSubrequestParams) -> np.ndarray:
        if hasattr(self, 'tokemak_overrides'):
            if params.tokemak in self.tokemak_overrides:
                return self.tokemak_overrides[params.tokemak](params)
        return self._get_signal_domain(params)
    
    @abstractmethod
    def _get_signal_domain(self, params : SignalDomainSubrequestParams) -> np.ndarray:
        """
        Default implementation of get_signal_domain.
        Used for any tokemak not in tokemak_overrides.
        """
        pass
    
class FlattopSubrequest(SignalDomainSubrequest):
    
    
# class SignalDomainSubrequest(Enum):
#     FULL = "full"
#     FLATTOP = "flattop"
#     RAMP_UP_AND_FLATTOP = "rampup_and_flattop"