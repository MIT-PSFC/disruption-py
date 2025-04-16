#!/usr/bin/env python3

"""
Module for generic physics methods.
"""

from typing import Dict

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak
from disruption_py.machine.cmod import CmodPhysicsMethods
from disruption_py.machine.d3d import D3DPhysicsMethods
from disruption_py.machine.east import EastPhysicsMethods


class GenericPhysicsMethods:
    """
    Class to hold generic physics methods.
    """

    @staticmethod
    @physics_method(columns=["shot_domain"])
    def get_shot_domain(params: PhysicsMethodParams):
        r"""
        Get the domain (or phase) of every time point in a shot and return it
        as a categorical feature:

        - 1: ramp-up
        - 0: flat-top
        - -1: ramp-down

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the categorical feature `shot_domain`.

        References
        -------
        - original source:
            - cmod: [get_flattop_indices.m](https://github.com/MIT-PSFC/disruption-py
            /blob/matlab/CMOD/matlab-core/get_flattop_indices.m)
            - d3d: [get_flattop_indices.m](https://github.com/MIT-PSFC/disr
            uption-py/blob/matlab/DIII-D/get_flattop_indices.m)
            - east: [get_flattop_indices.m](https:/github.com/MIT-PSFC/disruption-py/
            blob/matlab/EAST/utils/get_flattop_indices.m), [get_flattop_times.m](https://github
            .com/MIT-PSFC/disruption-py/blob/matlab/EAST/utils/get_flattop_times.m)
        - pull requests: #[433](https://github.com/MIT-PSFC/disruption-py/pull/433)
        - issues: #[408](https://github.com/MIT-PSFC/disruption-py/issues/408)
        """
        # Initialize dictionaries
        thresholds = {"dipprog_dt": None, "ip_prog": None, "power_supply_railed": None}
        signals = {"dipprog_dt": None, "ip_prog": None, "power_supply_railed": None}
        conditions = {
            "dipprog_dt": lambda signal, threshold: np.abs(signal) <= threshold,
            "ip_prog": lambda signal, threshold: np.abs(signal) >= threshold,
            "power_supply_railed": lambda signal, railed: signal != railed,
        }
        # Get data and threshold parameters
        if params.tokamak == Tokamak.CMOD:
            thresholds["dipprog_dt"] = 50e3
            thresholds["ip_prog"] = 100e3
            ip_parameters = CmodPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
            signals["ip_prog"] = ip_parameters["ip_prog"]
        elif params.tokamak == Tokamak.D3D:
            thresholds["dipprog_dt"] = 2e3
            thresholds["ip_prog"] = 100e3
            thresholds["power_supply_railed"] = 1
            ip_parameters = D3DPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
            signals["ip_prog"] = ip_parameters["ip_prog"]
        elif params.tokamak == Tokamak.EAST:
            thresholds["dipprog_dt"] = 1e3
            ip_parameters = EastPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
        else:
            return {"shot_domain": [np.nan]}

        shot_domain = np.full(len(signals["dipprog_dt"]), np.nan)
        # Get flattop domain indices
        indices_flattop = np.arange(len(shot_domain))
        for name in ["dipprog_dt", "ip_prog", "power_supply_railed"]:
            if all(v is not None for v in (signals[name], thresholds[name])):
                (indices,) = np.where(conditions[name](signals[name], thresholds[name]))
                indices_flattop = np.intersect1d(
                    indices_flattop, indices, assume_unique=True
                )

        # Get the longest subsequence of indices_flattop
        indices_flattop = max(
            np.split(indices_flattop, np.where(np.diff(indices_flattop) != 1)[0] + 1),
            key=len,
        )
        # Assign shot domains
        if len(indices_flattop) == 0:
            # Shot only has ramp up phase
            shot_domain[:] = 1
        else:
            flattop_start, flattop_end = indices_flattop[0], indices_flattop[-1] + 1
            shot_domain[:flattop_start] = 1
            shot_domain[flattop_start:flattop_end] = 0
            shot_domain[flattop_end:] = -1

        return {"shot_domain": shot_domain}
