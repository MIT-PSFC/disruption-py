#!/usr/bin/env python3

"""
Module initialization for MAST.
"""

from disruption_py.machine.mast.efit import MastEfitMethods
from disruption_py.machine.mast.physics import MastPhysicsMethods

METHOD_HOLDERS = [MastPhysicsMethods, MastEfitMethods]
