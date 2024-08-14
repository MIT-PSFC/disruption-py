#!/usr/bin/env python3

from disruption_py.machine.cmod.efit import CmodEfitMethods
from disruption_py.machine.cmod.physics import CmodPhysicsMethods
from disruption_py.machine.cmod.tearing import CmodTearingMethods

METHOD_HOLDERS = [CmodPhysicsMethods, CmodEfitMethods, CmodTearingMethods]
