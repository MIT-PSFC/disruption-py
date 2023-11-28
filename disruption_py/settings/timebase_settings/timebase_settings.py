from dataclasses import dataclass
from typing import Union, Callable
from disruption_py.settings.timebase_settings.set_times_requests import SetTimesRequest
from disruption_py.utils.mappings.mappings_helpers import map_string_attributes_to_enum
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
    
class SignalDomain(Enum):
    FULL = "full"
    FLATTOP = "flattop"
    RAMP_UP_AND_FLATTOP = "rampup_and_flattop"
    
@dataclass
class TimebaseSettings:
    set_times_request : SetTimesRequest = "efit"
    signal_domain : SignalDomain = "full"
    override_exising_data : bool = True
    interpolation_method : Union[str, Callable] = "linear"
    
    def __post_init__(self):
        map_string_attributes_to_enum(self, {
            "signal_type": SignalDomain,
            "interpolation_method": InterpolationMethod
        })
    