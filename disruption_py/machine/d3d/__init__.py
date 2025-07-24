#!/usr/bin/env python3

"""Package initialization for the DIII-D machine module."""

from disruption_py.machine.d3d.efit import D3DEfitMethods
from disruption_py.machine.d3d.physics import D3DPhysicsMethods
from disruption_py.machine.d3d.mirnov import D3DMirnovMethods

METHOD_HOLDERS = [D3DPhysicsMethods, D3DEfitMethods, D3DMirnovMethods]
