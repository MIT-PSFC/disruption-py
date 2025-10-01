#!/usr/bin/env python3

"""
Module for generic physics methods.
"""
import numpy as np

from disruption_py.config import config
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.cmod import CmodPhysicsMethods
from disruption_py.machine.d3d import D3DPhysicsMethods
from disruption_py.machine.east import EastPhysicsMethods
from disruption_py.machine.tokamak import Tokamak


class GenericPhysicsMethods:
    """
    Class to hold generic physics methods.
    """

    @staticmethod
    @physics_method(columns=["time_domain"])
    def get_time_domain(params: PhysicsMethodParams):
        r"""
        Get the domain (or phase) of every time point in a shot and return it
        as a categorical feature:

        - 1: ramp-up
        - 2: flat-top
        - 3: ramp-down

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the categorical feature `time_domain`.

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
        signals = {}
        thresholds = config(params.tokamak).physics.time_domain_thresholds
        conditions = {
            "dipprog_dt": lambda signal, threshold: np.abs(signal) <= threshold,
            "ip_prog": lambda signal, threshold: np.abs(signal) >= threshold,
            "power_supply_railed": lambda signal, railed: signal != railed,
        }
        # Get data and threshold parameters
        if params.tokamak == Tokamak.CMOD:
            ip_parameters = CmodPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
            signals["ip_prog"] = ip_parameters["ip_prog"]
        elif params.tokamak == Tokamak.D3D:
            ip_parameters = D3DPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
            signals["ip_prog"] = ip_parameters["ip_prog"]
            signals["power_supply_railed"] = ip_parameters["power_supply_railed"]
        elif params.tokamak == Tokamak.EAST:
            ip_parameters = EastPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
        else:
            return {"time_domain": [np.nan]}

        time_domain = np.full(len(params.times), np.nan)
        # Get flattop domain indices
        indices_flattop = np.arange(len(time_domain))
        for name in ["dipprog_dt", "ip_prog", "power_supply_railed"]:
            sig, thr = signals.get(name, None), thresholds.get(name, None)
            if all(v is not None for v in (sig, thr)):
                (indices,) = np.where(conditions[name](sig, thr))
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
            time_domain[:] = 1
        else:
            flattop_start, flattop_end = indices_flattop[0], indices_flattop[-1] + 1
            time_domain[:flattop_start] = 1
            time_domain[flattop_start:flattop_end] = 2
            time_domain[flattop_end:] = 3

        return {"time_domain": time_domain}
