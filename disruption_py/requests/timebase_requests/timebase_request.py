from dataclasses import dataclass
from typing import Union, Callable
from disruption_py.requests.timebase_requests.times_subrequests import SetTimesSubrequest
from disruption_py.utils.mappings.mappings_helpers import map_string_attributes_to_enum
from disruption_py.requests.timebase_requests.signal_domain_subrequest import SignalType
from enum import Enum

class InterpolationMethod(Enum):
    LINEAR = "linear"
    LOG = "log"
    LOG_LINEAR = "log_linear"
    EXP = "exp"
    EXP_LINEAR = "exp_linear"
    SIN = "sin"
    SIN_LINEAR = "sin_linear"
    SIN_EXP = "sin_exp"
    
@dataclass
class TimebaseRequest:
    times_request : SetTimesSubrequest
    signal_type : SignalType
    override_exising_data : bool = False
    interpolation_method : Union[str, Callable] = "linear"
    
    def __post_init__(self):
        map_string_attributes_to_enum(self, {
            "signal_type": SignalType,
            "interpolation_method": InterpolationMethod
        })
    