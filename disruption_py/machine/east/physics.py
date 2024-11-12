#!/usr/bin/env python3

"""
Module for retrieving and calculating data for DIII-D physics methods.
"""

import traceback

import numpy as np
import scipy
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import (
    interp1,
    matlab_get_bolo,
    matlab_gsastd,
    matlab_power,
)
from disruption_py.machine.tokamak import Tokamak


class EASTPhysicsMethods:
    """
    A class to retrieve and calculate physics-related data for DIII-D.
    """
