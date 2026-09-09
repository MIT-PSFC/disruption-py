#!/usr/bin/env python3

"""
Module for generic physics methods.
"""

import numpy as np

from disruption_py.config import config
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.cmod import CmodPhysicsMethods
from disruption_py.machine.d3d import D3DPhysicsMethods
from disruption_py.machine.east import EastPhysicsMethods
from disruption_py.machine.east.util import EastUtilMethods
from disruption_py.machine.generic.util import GenericUtilMethods
from disruption_py.machine.mast.physics import MastPhysicsMethods
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
            The parameters containing the data connection, shot id and more.

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
            psr_parameter = D3DPhysicsMethods.get_power_supply_railed(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
            signals["ip_prog"] = ip_parameters["ip_prog"]
            signals["power_supply_railed"] = psr_parameter["power_supply_railed"]
        elif params.tokamak == Tokamak.EAST:
            ip_parameters = EastPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
        elif params.tokamak == Tokamak.MAST:
            ip_parameters = MastPhysicsMethods.get_ip_parameters(params=params)
            signals["dipprog_dt"] = ip_parameters["dipprog_dt"]
            signals["ip_prog"] = ip_parameters["ip_prog"]
        else:
            return {"time_domain": [np.nan]}

        # Check if all signals are available and valid
        for signal in signals.values():
            if np.isnan(signal).all():
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

    @staticmethod
    @physics_method(columns=["current_quench_time"])
    def get_current_quench_time(params: PhysicsMethodParams):
        """
        Determine and compute the current quench time of a shot. If a shot is determined
        to be non-disruptive, return the current quench time as NaN. The criteria for a
        disruptive shot are as follows:

        - 1. Shot duration > duration_min: reject very short shots.
        - 2. abs(Ip_max) > ip_threshold: reject very low current shots.
        - 3. Ip0 / Ip_max > rampdown_threshold: reject shots that disrupt late in current ramp down.
        - 4. -Ip0 / max_dIdt < tau_CQ_max: reject shots with relatively slow current decay.
        - 5. abs(Ip_final) < Ip_final_max: reject minor disruptions.
        - 6. abs(Ip0) > Ip_threshold: reject very low current disruptions.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the data connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the `current_quench_time`.

        References
        -------
        - original source:
            - cmod: [test_disrupt.pro]
            - d3d: [test_for_disruption.m]
            - east: [test_for_disruption.m]
            - kstar: [test_for_disruption_kstar.m]
        - pull requests: #[545](https://github.com/MIT-PSFC/disruption-py/pull/545)
        - issues: #[223](https://github.com/MIT-PSFC/disruption-py/issues/223)
        """
        # Initialize test criteria
        criteria = {
            "shot_duration": lambda duration, duration_min: duration > duration_min,
            "abs_ip_max": lambda ip_max, ip_threshold: abs(ip_max) > ip_threshold,
            "ip0_over_ip_max": lambda ip0_over_ip_max, rampdown_threshold: ip0_over_ip_max
            > rampdown_threshold,
            "ip0_over_max_didt": lambda ip0_over_max_didt, tau_cq_max: -ip0_over_max_didt
            < tau_cq_max,
            "abs_ip_final": lambda ip_final, ip_final_max: abs(ip_final) < ip_final_max,
            "abs_ip0": lambda ip0, ip_threshold: abs(ip0) > ip_threshold,
        }

        # Get ip and ip timebase
        if params.tokamak == Tokamak.CMOD:
            ip, t_ip = params.get_data_with_dims(
                r"\ip", tree_name="magnetics"
            )  # [A], [s]
        elif params.tokamak == Tokamak.D3D:
            ip, t_ip = params.get_data_with_dims(
                f"ptdata('ip', {params.shot_id})"
            )  # [A], [ms]
            t_ip = t_ip / 1.0e3  # [ms] -> [s]
            # Subtract baseline offset from ip
            (baseline_indices,) = np.where(t_ip <= 0)
            if len(baseline_indices) > 0:
                ip_baseline = np.mean(ip[baseline_indices])
                ip -= ip_baseline
        elif params.tokamak == Tokamak.EAST:
            ip, t_ip = EastUtilMethods.retrieve_ip(params, params.shot_id)
        elif params.tokamak == Tokamak.HBTEP:
            ip, t_ip = params.get_data_with_dims(
                r"\top.sensors.rogowskis:ip", tree_name="hbtep2"
            )  # [A], [s]
        elif params.tokamak == Tokamak.MAST:
            ip = params.get_data("summary/ip")
            t_ip = params.get_data("summary/time")
        else:
            raise NotImplementedError

        # Get test thresholds
        thresholds = config(params.tokamak).physics.current_quench_time_thresholds

        # Compute parameters for the tests
        # Get the duration and polarity
        end_of_current_params = GenericUtilMethods.get_end_of_current(
            ip=ip, ip_time=t_ip, threshold=thresholds["abs_ip_max"]
        )
        duration = end_of_current_params["duration"]
        polarity = end_of_current_params["polarity"]
        # Check for incomplete digitization
        if duration < 0:
            raise CalculationError(
                "Plasma current signal may not have been digitized for long enough "
                "to capture the end of the discharge. Cannot determine current_quench_time."
            )
        # Check for a bona-fide plasma
        if duration == 0:
            params.logger.verbose(
                "get_current_quench_time: Shot does not satisfy the "
                "bona fide plasma criterion. Marking as non-disruptive."
            )
            return {"current_quench_time": [np.nan]}

        # Find the maximum plasma current excluding the current spike
        ip_upright = ip * polarity
        (time_indices,) = np.where((t_ip > 0) & (t_ip < duration - 0.050))
        if len(time_indices) == 0:
            raise CalculationError(
                "Invalid timebase. Cannot compute current quench time."
            )
        ip_max = max(ip_upright[time_indices]) * polarity
        if ip_max == 0:
            # ip_max is in the denominator of `ip0_over_ip_max` so it cannot be 0.
            raise CalculationError(
                "Maximum plasma current is zero. Cannot compute current quench time."
            )

        # Find the plasma current right before the end of the discharge
        (time_indices,) = np.where((t_ip > duration - 0.06) & (t_ip < duration - 0.04))
        if len(time_indices) == 0:
            raise CalculationError(
                "Invalid timebase. Cannot compute current quench time."
            )
        candidate_ip0 = np.mean(ip_upright[time_indices]) * polarity

        # Compute dI/dt during the latter part of the discharge
        (time_indices,) = np.where((t_ip > duration - 0.05) & (t_ip < duration + 0.05))
        didt_upright = np.diff(ip_upright[time_indices]) / np.diff(t_ip[time_indices])
        if len(didt_upright) == 0:
            raise CalculationError(
                "Invalid timebase. Cannot compute current quench time."
            )
        indx = np.argmin(didt_upright)
        candidate_max_didt = didt_upright[indx] * polarity
        candidate_t_disrupt = t_ip[time_indices[indx]]
        if candidate_max_didt == 0:
            # candidate_max_didt is in the denominator of `ip0_over_max_didt` so it cannot be 0.
            raise CalculationError(
                "Maximum dI/dt is zero. Cannot compute current quench time."
            )

        # Compute ip_final
        (time_indices,) = np.where(
            (t_ip > candidate_t_disrupt) & (t_ip < candidate_t_disrupt + 0.15)
        )
        if len(time_indices) == 0:
            raise CalculationError(
                "Invalid timebase. Cannot compute current quench time."
            )
        ip_final = min(ip_upright[time_indices])

        # Collect the computed parameters for the tests
        parameters = {
            "shot_duration": duration,
            "abs_ip_max": ip_max,
            "ip0_over_ip_max": candidate_ip0 / ip_max,
            "ip0_over_max_didt": candidate_ip0 / candidate_max_didt,
            "abs_ip_final": ip_final,
            "abs_ip0": candidate_ip0,
        }

        # Run all 6 tests
        for test, criterion in criteria.items():
            parameter, threshold = parameters[test], thresholds[test]
            # Use threshold = None to indicate that we will skip a test.
            if threshold is None:
                continue
            # Run the test criteriion.
            if not criterion(parameter, threshold):
                # Failed the test, mark the shot as non-disruptive.
                return {"current_quench_time": [np.nan]}

        return {"current_quench_time": np.full(len(params.times), candidate_t_disrupt)}
