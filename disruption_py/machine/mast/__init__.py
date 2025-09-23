#!/usr/bin/env python3

"""
Module initialization for MAST.
"""

from disruption_py.machine.mast.physics import MastPhysicsMethods
from disruption_py.machine.mast.efit import MastEfitMethods

METHOD_HOLDERS = [MastPhysicsMethods, MastEfitMethods]
