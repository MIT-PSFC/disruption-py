#!/usr/bin/env python3

from disruption_py.machine.cmod.efit import CmodEfitMethods
from disruption_py.machine.cmod.physics import CmodPhysicsMethods, CmodTearingMethods

METHOD_HOLDERS = [CmodPhysicsMethods, CmodEfitMethods, CmodTearingMethods]
