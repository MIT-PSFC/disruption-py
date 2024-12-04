#!/usr/bin/env python3

"""Package initialization for the EAST machine module."""

from disruption_py.machine.east.efit import EastEfitMethods

from disruption_py.machine.east.physics import EastPhysicsMethods

METHOD_HOLDERS = [EastEfitMethods, EastPhysicsMethods]
