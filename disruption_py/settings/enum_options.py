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