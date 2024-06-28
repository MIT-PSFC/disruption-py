#!/usr/bin/env python3

from disruption_py.shots.parameter_methods.d3d.basic_parameter_methods import (
    BasicD3DRequests,
)
from disruption_py.shots.parameter_methods.d3d.efit_parameter_methods import (
    D3DEfitRequests,
)

D3D_DEFAULT_PARAMETER_METHODS = [
    D3DEfitRequests(),
    BasicD3DRequests(),
]
