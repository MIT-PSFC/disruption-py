#!/usr/bin/env python3

from disruption_py.machine.cmod.basic import (
    BasicCmodRequests,
    CModEfitRequests,
)

CMOD_DEFAULT_PARAMETER_METHODS = [
    CModEfitRequests(),
    BasicCmodRequests(),
]
