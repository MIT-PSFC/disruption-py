#!/usr/bin/env python3

"""Package initialization for the CMOD machine module."""

from disruption_py.machine.cmod.efit import CmodEfitMethods
from disruption_py.machine.cmod.physics import CmodPhysicsMethods
from disruption_py.machine.cmod.mirnov import CmodMirnovMethods

METHOD_HOLDERS = [CmodPhysicsMethods, CmodEfitMethods, CmodMirnovMethods]
