#!/usr/bin/env python3

from disruption_py.shots.parameter_methods.cmod.basic_parameter_methods import (
    BasicCmodRequests,
    CModEfitRequests,
)

CMOD_DEFAULT_PARAMETER_METHODS = [
    CModEfitRequests(),
    BasicCmodRequests(),
]
