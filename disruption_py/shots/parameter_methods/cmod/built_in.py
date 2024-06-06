#!/usr/bin/env python3

from disruption_py.shots.parameter_methods.cmod.basic_parameter_methods import (
    BasicCmodRequests,
    CModEfitRequests,
)

CMOD_DEFAULT_SHOT_DATA_REQUESTS = [
    CModEfitRequests(),
    BasicCmodRequests(),
]
