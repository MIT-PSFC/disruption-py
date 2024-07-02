#!/usr/bin/env python3

from disruption_py.machine.d3d.basic import (
    BasicD3DRequests,
)
from disruption_py.machine.d3d.efit import (
    D3DEfitRequests,
)

D3D_DEFAULT_PARAMETER_METHODS = [
    D3DEfitRequests(),
    BasicD3DRequests(),
]
