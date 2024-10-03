#!/usr/bin/env python3

"""
Module for retrieving and calculating data for CMOD physics methods.
"""

import traceback
import warnings

import numpy as np
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import (
    gaussian_fit,
    gaussian_fit_with_fixed_mean,
    interp1,
    smooth,
)
from disruption_py.machine.cmod.thomson import CmodThomsonDensityMeasure
from disruption_py.machine.tokamak import Tokamak

warnings.filterwarnings("error", category=RuntimeWarning)


class CmodPhysicsMethods:
    """
    This class provides methods to retrieve and calculate physics-related data
    for CMOD.
    """

    @staticmethod
    @cache_method
    def _get_active_wire_segments(params: PhysicsMethodParams):
        """
        Retrieve active wire segments from the MDSplus tree.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection and shot info.

        Returns
        -------
        list of tuple
            A list of tuples, where each tuple contains the node path of the
            active segment and its start time. The list is sorted by start time.
        """
        params.mds_conn.open_tree(tree_name="pcs")
        root_nid = params.mds_conn.get("GetDefaultNid()")
        children_nids = params.mds_conn.get(
            'getnci(getnci($, "CHILDREN_NIDS"), "NID_NUMBER")', arguments=root_nid
        )
        children_paths = params.mds_conn.get(
            'getnci($, "FULLPATH")', arguments=children_nids
        )

        # Collect active segments and their information
        active_segments = []
        for node_path in children_paths:
            node_path = node_path.strip()
            if node_path.split(".")[-1].startswith("SEG_"):
                is_on = params.mds_conn.get_data(
                    'getnci($, "STATE")', arguments=node_path + ":SEG_NUM"
                )
                # 0 represents node being on, 1 represents node being off
                if is_on != 0:
                    continue
                active_segments.append(
                    (
                        node_path,
                        params.mds_conn.get_data(
                            node_path + ":start_time", tree_name="pcs"
                        ),
                    )
                )

        active_segments.sort(key=lambda n: n[1])
        return active_segments

    @staticmethod
    @physics_method(columns=["time_until_disrupt"], tokamak=Tokamak.CMOD)
    def get_time_until_disrupt(params: PhysicsMethodParams):
        """
        Calculate the time until disruption.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the disruption information and times.

        Returns
        -------
        dict
            A dictionary with a single key "time_until_disrupt" containing a list
            of time until disruption.
        """
        time_until_disrupt = [np.nan]
        if params.disrupted:
            time_until_disrupt = params.disruption_time - params.times
        return {"time_until_disrupt": time_until_disrupt}

    @staticmethod
    def _get_ip_parameters(times, ip, magtime, ip_prog, pcstime):
        """
        Calculates actual and programmed current as well as their derivatives
        and difference.

        The time derivatives are useful for discriminating between rampup, flattop,
        and rampdown.

        Parameters
        ----------
        times : array_like
            Time array for the shot.
        ip : array_like
            Actual plasma current.
        magtime : array_like
            Time array for the plasma current.
        ip_prog : array_like
            Programmed plasma current.
        pcstime : array_like
            Time array for the programmed plasma current.

        Returns
        -------
        ip : array_like
            Actual plasma current.
        dip_dt : array_like
            Time derivative of the actual plasma current.
        dip_smoothed : array_like
            Smoothed time derivative of the actual plasma current.
        ip_prog : array_like
            Programmed plasma current.
        dipprog_dt : array_like
            Time derivative of the programmed plasma current.
        ip_error : array_like
            Difference between the actual and programmed plasma current.

        Original Authors
        ----------------
        - Alex Tinguely
        - Robert Granetz
        - Ryan Sweeney

        Sources
        -------
        - matlab/cmod_matlab/matlab-core/get_Ip_parameters.m
        - matlab/cmod_matlab/matlab-core/get_Ip_parameters.m
        """
        dip = np.gradient(ip, magtime)
        dip_smoothed = smooth(dip, 11)  # ,ends_type=0)
        dipprog_dt = np.gradient(ip_prog, pcstime)
        ip_prog = interp1(
            pcstime, ip_prog, times, bounds_error=False, fill_value=ip_prog[-1]
        )
        dipprog_dt = interp1(pcstime, dipprog_dt, times, bounds_error=False)
        ip = interp1(magtime, ip, times)
        dip = interp1(magtime, dip, times)
        dip_smoothed = interp1(magtime, dip_smoothed, times)

        ip_error = (np.abs(ip) - np.abs(ip_prog)) * np.sign(ip)
        # import pdb; pdb.set_trace()
        output = {
            "ip": ip,
            "dip_dt": dip,
            "dip_smoothed": dip_smoothed,
            "ip_prog": ip_prog,
            "dipprog_dt": dipprog_dt,
            "ip_error": ip_error,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["ip", "dip_dt", "dip_smoothed", "ip_prog", "dipprog_dt", "ip_error"],
        tokamak=Tokamak.CMOD,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """
        Retrieve and interpolate Ip parameters.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the interpolated Ip parameters, including
            "ip", "dip_dt", "dip_smoothed", "ip_prog", "dipprog_dt", and "ip_error".
        """
        # Automatically generated
        active_segments = CmodPhysicsMethods._get_active_wire_segments(params=params)

        # Default PCS timebase is 1 KHZ
        pcstime = np.array(np.arange(-4, 12.383, 0.001))
        ip_prog = np.full(pcstime.shape, np.nan)

        # For each active segment:
        # 1.) Find the wire for IP control and check if it has non-zero PID gains
        # 2.) IF it does, interpolate IP programming onto the PCS timebase
        # 3.) Clip to the start and stop times of PCS timebase
        for node_path, start in active_segments:
            # Ip wire can be one of 16 but is normally no. 16
            for wire_index in range(16, 0, -1):
                wire_node_name = params.mds_conn.get_data(
                    node_path + f":P_{wire_index :02d}:name", tree_name="pcs"
                )
                if wire_node_name == "IP":
                    try:
                        pid_gains = params.mds_conn.get_data(
                            node_path + f":P_{wire_index :02d}:pid_gains",
                            tree_name="pcs",
                        )
                        if np.any(pid_gains):
                            signal, sigtime = params.mds_conn.get_data_with_dims(
                                node_path + f":P_{wire_index :02d}", tree_name="pcs"
                            )
                            ip_prog_temp = interp1(
                                sigtime,
                                signal,
                                pcstime,
                                bounds_error=False,
                                fill_value=signal[-1],
                            )
                            end = pcstime[
                                np.argmin(np.abs(pcstime - sigtime[-1]) + 0.0001)
                            ]
                            segment_indices = np.where(
                                (pcstime >= start) & (pcstime <= end)
                            )
                            ip_prog[segment_indices] = ip_prog_temp[segment_indices]
                    except mdsExceptions.MdsException:
                        params.logger.warning(
                            "[Shot %s]: Error getting PID gains for wire %s",
                            params.shot_id,
                            wire_index,
                        )
                        params.logger.debug(
                            "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                        )
                    break  # Break out of wire_index loop
        ip, magtime = params.mds_conn.get_data_with_dims(
            r"\ip", tree_name="magnetics", astype="float64"
        )
        output = CmodPhysicsMethods._get_ip_parameters(
            params.times, ip, magtime, ip_prog, pcstime
        )
        return output

    @staticmethod
    def _get_z_parameters(times, z_prog, pcstime, z_error_without_ip, ip, dpcstime):
        """
        Get values of Z_error, Z_prog, and derived signals from plasma control
        system (PCS).

        Z_prog is the programmed vertical position of the plasma current centroid,
        and Z_error is the difference between the actual position and that requested
        (Z_error = Z_cur - Z_prog). Thus, the actual (estimated) position, Z_cur,
        can be calculated. And the vertical velocity, v_z, can be taken from the
        time derivative, and the product z_times_v_z ( = Z_cur * v_z) is also calculated.

        Parameters
        ----------
        times : array_like
            Time array for the shot.
        z_prog : array_like
            Programmed vertical position of the plasma current centroid.
        pcstime : array_like
            Time array for the programmed vertical position of the plasma current
            centroid.
        z_error_without_ip : array_like
            Difference between the actual and programmed vertical position of the
            plasma current centroid.
        ip : array_like
            Actual plasma current.
        dpcstime : array_like
            Time array for the actual plasma current.

        Returns
        -------
        z_error : array_like
            Difference between the actual and programmed vertical position of the
            plasma current centroid.
        z_prog : array_like
            Programmed vertical position of the plasma current centroid.
        z_cur : array_like
            Actual (estimated) vertical position of the plasma current centroid.
        v_z : array_like
            Vertical velocity.
        z_times_v_z : array_like
            Product of the vertical position and vertical velocity.

        Original Authors
        ----------------
        - Alex Tinguely
        - Robert Granetz

        Sources
        -------
        - matlab/cmod_matlab/matlab-core/get_Z_parameters.m
        - matlab/cmod_matlab/matlab-core/get_Z_parameters.m

        """
        divsafe_ip = np.where(ip != 0, ip, np.nan)
        z_error = -1 * z_error_without_ip / divsafe_ip  # [m]
        z_prog_dpcs = interp1(pcstime, z_prog, dpcstime)
        z_cur = z_prog_dpcs + z_error  # [m]
        v_z = np.gradient(z_cur, dpcstime)  # m/s
        z_times_v_z = z_cur * v_z  # m^2/s
        z_prog = interp1(pcstime, z_prog, times, "linear", False, z_prog[-1])
        z_error = interp1(dpcstime, z_error, times, "linear", False, z_error[-1])
        z_cur = interp1(dpcstime, z_cur, times, "linear", False, z_cur[-1])
        v_z = interp1(dpcstime, v_z, times, "linear", False, v_z[-1])
        z_times_v_z = interp1(
            dpcstime, z_times_v_z, times, "linear", False, z_times_v_z[-1]
        )
        output = {
            "z_error": z_error,
            "z_prog": z_prog,
            "zcur": z_cur,
            "v_z": v_z,
            "z_times_v_z": z_times_v_z,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["z_error", "z_prog", "zcur", "v_z", "z_times_v_z"],
        tokamak=Tokamak.CMOD,
    )
    def get_z_parameters(params: PhysicsMethodParams):
        """
        Retrieve and interpolate plasma's vertical position parameters.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the vertical position parameters, including "z_error", "z_prog",
            "zcur", "v_z", and "z_times_v_z".
        """
        pcstime = np.array(np.arange(-4, 12.383, 0.001))
        z_prog = np.empty(pcstime.shape)
        z_prog.fill(np.nan)
        z_prog_temp = z_prog.copy()
        z_wire_index = -1
        active_wire_segments = CmodPhysicsMethods._get_active_wire_segments(
            params=params
        )

        for node_path, start in active_wire_segments:
            for wire_index in range(1, 17):
                wire_node_name = params.mds_conn.get_data(
                    node_path + f":P_{wire_index :02d}:name", tree_name="pcs"
                )
                if wire_node_name == "ZCUR":
                    try:
                        pid_gains = params.mds_conn.get_data(
                            node_path + f":P_{wire_index :02d}:pid_gains",
                            tree_name="pcs",
                        )
                        if np.any(pid_gains):
                            signal, sigtime = params.mds_conn.get_data_with_dims(
                                node_path + f":P_{wire_index :02d}", tree_name="pcs"
                            )
                            end = sigtime[
                                np.argmin(np.abs(sigtime - pcstime[-1]) + 0.0001)
                            ]
                            z_prog_temp = interp1(
                                sigtime,
                                signal,
                                pcstime,
                                "linear",
                                False,
                                fill_value=signal[-1],
                            )
                            z_wire_index = wire_index
                            segment_indices = [
                                np.where((pcstime >= start) & (pcstime <= end))
                            ]
                            z_prog[segment_indices] = z_prog_temp[segment_indices]
                            break
                    except mdsExceptions.MdsException:
                        params.logger.debug(
                            "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                        )
                        continue  # TODO: Consider raising appropriate error
                else:
                    continue
                break
        if z_wire_index == -1:
            raise CalculationError("Data source error: No ZCUR wire was found")
        # Read in A_OUT, which is a 16xN matrix of the errors for *all* 16 wires for
        # *all* of the segments. Note that DPCS time is usually taken at 10kHz.
        wire_errors, dpcstime = params.mds_conn.get_data_with_dims(
            r"\top.hardware.dpcs.signals:a_out", tree_name="hybrid", dim_nums=[1]
        )
        # The value of Z_error we read is not in the units we want. It must be *divided*
        #  by a factor AND *divided* by the plasma current.
        z_error_without_factor_and_ip = wire_errors[:, z_wire_index]
        z_error_without_ip = np.empty(z_error_without_factor_and_ip.shape)
        z_error_without_ip.fill(np.nan)
        # Also, it turns out that different segments have different factors. So we
        # search through the active segments (determined above), find the factors,
        # and *divide* by the factor only for the times in the active segment (as
        # determined from start_times and stop_times.
        for i, (_, start) in enumerate(active_wire_segments):
            if i == len(active_wire_segments) - 1:
                end = pcstime[-1]
            else:
                end = active_wire_segments[i + 1][1]
            z_factor = params.mds_conn.get_data(
                rf"\dpcs::top.seg_{i+1:02d}:p_{z_wire_index:02d}:predictor:factor",
                tree_name="hybrid",
            )
            temp_indx = np.where((dpcstime >= start) & (dpcstime <= end))
            z_error_without_ip[temp_indx] = (
                z_error_without_factor_and_ip[temp_indx] / z_factor
            )  # [A*m]
        # Next we grab ip, which comes from a_in:input_056. This also requires
        # *multiplication* by a factor.
        # NOTE that I can't get the following ip_without_factor to work for shots
        # before 2015.
        # TODO: Try to fix this
        if params.shot_id > 1150101000:
            ip_without_factor = params.mds_conn.get_data(
                r"\hybrid::top.hardware.dpcs.signals.a_in:input_056", tree_name="hybrid"
            )
            ip_factor = params.mds_conn.get_data(
                r"\hybrid::top.dpcs_config.inputs:input_056:p_to_v_expr",
                tree_name="hybrid",
            )
            ip = ip_without_factor * ip_factor  # [A]
        else:
            ip, ip_time = params.mds_conn.get_data_with_dims(
                r"\ip", tree_name="magnetics"
            )
            ip = interp1(ip_time, ip, dpcstime)
        return CmodPhysicsMethods._get_z_parameters(
            params.times, z_prog, pcstime, z_error_without_ip, ip, dpcstime
        )

    @staticmethod
    def _get_ohmic_parameters(
        times, v_loop, v_loop_time, li, efittime, dip_smoothed, ip, r0
    ):
        """
        Calculate the ohmic power from the loop voltage, inductive voltage, and
        plasma current.

        Parameters
        ----------
        times : array_like
            The times at which to calculate the ohmic power.
        v_loop : array_like
            The loop voltage.
        v_loop_time : array_like
            The times at which the loop voltage was measured.
        li : array_like
            The plasma's internal inductance from EFIT.
        efittime : array_like
            The EFIT time base.
        dip_smoothed : array_like
            The smoothed time derivative of the measured plasma current.
        ip : array_like
            The plasma current.
        r0 : array_like
            The major radius of the plasma's magnetic axis.

        Returns
        -------
        p_oh : array_like
            The ohmic power.
        v_loop : array_like
            The loop voltage.
        """
        inductance = 4.0 * np.pi * 1.0e-7 * r0 * li / 2.0
        v_loop = interp1(v_loop_time, v_loop, times)
        inductance = interp1(efittime, inductance, times)
        v_inductive = inductance * dip_smoothed
        v_resistive = v_loop - v_inductive
        p_ohm = ip * v_resistive
        output = {"p_oh": p_ohm, "v_loop": v_loop}
        return output

    @staticmethod
    @physics_method(
        columns=["p_oh", "v_loop"],
        tokamak=Tokamak.CMOD,
    )
    def get_ohmic_parameters(params: PhysicsMethodParams):
        """
        Retrieve and calculate ohmic heating parameters.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the calculated ohmic parameters, including
            "p_oh" and "v_loop".
        """
        v_loop, v_loop_time = params.mds_conn.get_data_with_dims(
            r"\top.mflux:v0", tree_name="analysis", astype="float64"
        )
        if len(v_loop_time) <= 1:
            raise CalculationError("No data for v_loop_time")

        li, efittime = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:li", tree_name="_efit_tree", astype="float64"
        )  # [dimensionless], [s]
        ip_parameters = CmodPhysicsMethods.get_ip_parameters(params=params)
        r0 = 0.01 * params.mds_conn.get_data(
            r"\efit_aeqdsk:rmagx", tree_name="_efit_tree"
        )  # [cm] -> [m]

        output = CmodPhysicsMethods._get_ohmic_parameters(
            params.times,
            v_loop,
            v_loop_time,
            li,
            efittime,
            ip_parameters["dip_smoothed"],
            ip_parameters["ip"],
            r0,
        )
        return output

    @staticmethod
    def _get_power(times, p_lh, t_lh, p_icrf, t_icrf, p_rad, t_rad, p_ohm):
        """
        Calculate the input and radiated powers, and then calculate the
        radiated fraction.

        Parameters
        ----------
        times : np.ndarray
            The time array for which to calculate the power.
        p_lh : np.ndarray or None
            The power from lower hybrid heating.
        t_lh : np.ndarray or None
            The time array corresponding to lower hybrid heating power.
        p_icrf : np.ndarray or None
            The power from ICRF heating.
        t_icrf : np.ndarray or None
            The time array corresponding to ICRF heating power.
        p_rad : np.ndarray or None
            The radiated power.
        t_rad : np.ndarray or None
            The time array corresponding to radiated power.
        p_ohm : np.ndarray
            The ohmic heating power.

        Returns
        -------
        dict
            A dictionary containing the calculated power values, including
            "p_rad", "dprad_dt", "p_lh", "p_icrf", "p_input", and "radiated_fraction".
        """
        if p_lh is not None and isinstance(t_lh, np.ndarray) and len(t_lh) > 1:
            p_lh = interp1(t_lh, p_lh * 1.0e3, times)
        else:
            p_lh = np.zeros(len(times))

        if p_icrf is not None and isinstance(t_icrf, np.ndarray) and len(t_icrf) > 1:
            p_icrf = interp1(t_icrf, p_icrf * 1.0e6, times, bounds_error=False)
        else:
            p_icrf = np.zeros(len(times))

        if (
            t_rad is None
            or p_rad is None
            or not isinstance(t_rad, np.ndarray)
            or len(t_rad) <= 1
        ):
            p_rad = np.array([np.nan] * len(times))  # TODO: Fix
            dprad = p_rad.copy()
        else:
            p_rad = p_rad * 1.0e3  # [W]
            p_rad = p_rad * 4.5  # Factor of 4.5 comes from cross-calibration with
            # 2pi_foil during flattop times of non-disruptive
            # shots, excluding times for
            # which p_rad (uncalibrated) <= 1.e5 W
            dprad = np.gradient(p_rad, t_rad)
            p_rad = interp1(t_rad, p_rad, times)
            dprad = interp1(t_rad, dprad, times)
        p_input = p_ohm + p_lh + p_icrf
        rad_fraction = p_rad / p_input
        rad_fraction[rad_fraction == np.inf] = np.nan
        output = {
            "p_rad": p_rad,
            "dprad_dt": dprad,
            "p_lh": p_lh,
            "p_icrf": p_icrf,
            "p_input": p_input,
            "radiated_fraction": rad_fraction,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["p_rad", "dprad_dt", "p_lh", "p_icrf", "p_input", "radiated_fraction"],
        tokamak=Tokamak.CMOD,
    )
    def get_power(params: PhysicsMethodParams):
        """
        NOTE: the timebase for the LH power signal does not extend over the full
            time span of the discharge.  Therefore, when interpolating the LH power
            signal onto the "timebase" array, the LH signal has to be extrapolated
            with zero values.  This is an option in the 'interp1' routine.  If the
            extrapolation is not done, then the 'interp1' routine will assign NaN
            (Not-a-Number) values for times outside the LH timebase, and the NaN's
            will propagate into p_input and rad_fraction, which is not desirable.
        """
        values = [
            None
        ] * 6  # List to store the time and values of the LH power, icrf power, and radiated power
        trees = ["LH", "RF", "spectroscopy"]
        nodes = [r"\LH::TOP.RESULTS:NETPOW", r"\rf::rf_power_net", r"\twopi_diode"]
        for i in range(3):
            try:
                sig, sig_time = params.mds_conn.get_data_with_dims(
                    nodes[i], tree_name=trees[i], astype="float64"
                )
                values[2 * i] = sig
                values[2 * i + 1] = sig_time
            except (mdsExceptions.TreeFOPENR, mdsExceptions.TreeNNF):
                continue
        p_oh = CmodPhysicsMethods.get_ohmic_parameters(params=params)["p_oh"]
        output = CmodPhysicsMethods._get_power(params.times, *values, p_oh)
        return output

    @staticmethod
    def _get_kappa_area(times, aminor, area, a_times):
        """
        Calculate and interpolate kappa_area

        Parameters
        ----------
        times : np.ndarray
            The time array for which to calculate the kappa_area.
        aminor : np.ndarray
            The minor radius values.
        area : np.ndarray
            The area values.
        a_times : np.ndarray
            The time array corresponding to the area values.

        Returns
        -------
        dict
            A dictionary containing the kappa_area.
        """
        output = {"kappa_area": interp1(a_times, area / (np.pi * aminor**2), times)}
        return output

    @staticmethod
    @physics_method(columns=["kappa_area"], tokamak=Tokamak.CMOD)
    def get_kappa_area(params: PhysicsMethodParams):
        """
        Retrieve and calculate the plasma's ellipticity (kappa, also known as
        the elongation) using its area and minor radius.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the calculated "kappa_area".
        """
        aminor = params.mds_conn.get_data(
            r"\efit_aeqdsk:aminor", tree_name="_efit_tree", astype="float64"
        )
        area = params.mds_conn.get_data(
            r"\efit_aeqdsk:area", tree_name="_efit_tree", astype="float64"
        )
        times = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
        )

        aminor[aminor <= 0] = 0.001  # make sure aminor is not 0 or less than 0
        # make sure area is not 0 or less than 0
        area[area <= 0] = 3.14 * 0.001**2
        output = CmodPhysicsMethods._get_kappa_area(params.times, aminor, area, times)
        return output

    @staticmethod
    def _get_rotation_velocity(times, intensity, time, vel, hirextime):
        """
        Uses spectroscopy graphs of ionized (to hydrogen and helium levels) Argon
        to calculate velocity. Because of the heat profile of the plasma, suitable
        measurements are only found near the center.
        """
        v_0 = np.empty(len(time))
        # Check that the argon intensity pulse has a minimum count and duration
        # threshold
        valid_indices = np.where(intensity > 1000 & intensity < 10000)
        # Matlab code just multiplies by time delta but that doesn't work in the
        # case where we have different time deltas. Instead we sum the time deltas
        # for all valid indices to check the total duration
        if np.sum(time[valid_indices + 1] - time[valid_indices]) >= 0.2:
            v_0 = interp1(hirextime, vel, time)
            # TODO: Determine better threshold
            v_0[np.where(abs(v_0) > 200)] = np.nan
            v_0 *= 1000.0
        v_0 = interp1(time, v_0, times)
        return {"v_0": v_0}

    @staticmethod
    @physics_method(
        columns=["n_equal_1_mode", "n_equal_1_normalized", "n_equal_1_phase", "bt"],
        tokamak=Tokamak.CMOD,
    )
    def get_n_equal_1_amplitude(params: PhysicsMethodParams):
        """
        Calculate n=1 amplitude and phase.

        This method uses the four BP13 Bp sensors near the midplane on the outboard
        vessel wall.  The calculation is done by using a least squares fit to an
        expansion in terms of n = 0 & 1 toroidal harmonics.  The BP13 sensors are
        part of the set used for plasma control and equilibrium reconstruction,
        and their signals have been analog integrated (units: tesla), so they
        don't have to be numerically integrated.  These four sensors were working
        well in 2014, 2015, and 2016.  I looked at our locked mode MGI run on
        1150605, and the different applied A-coil phasings do indeed show up on
        the n=1 signal.

        N=1 toroidal assymmetry in the magnetic fields
        """
        # These sensors are placed toroidally around the machine. Letters refer to
        # the 2 ports the sensors were placed between.
        bp13_names = ["BP13_BC", "BP13_DE", "BP13_GH", "BP13_JK"]
        bp13_signals = np.empty((len(params.times), len(bp13_names)))

        path = r"\mag_bp_coils."
        bp_node_names = params.mds_conn.get_data(
            path + "nodename", tree_name="magnetics"
        )
        phi = params.mds_conn.get_data(path + "phi", tree_name="magnetics")
        btor_pickup_coeffs = params.mds_conn.get_data(
            path + "btor_pickup", tree_name="magnetics"
        )
        _, bp13_indices, _ = np.intersect1d(
            bp_node_names, bp13_names, return_indices=True
        )
        bp13_phi = phi[bp13_indices] + 360  # INFO
        bp13_btor_pickup_coeffs = btor_pickup_coeffs[bp13_indices]
        btor, t_mag = params.mds_conn.get_data_with_dims(
            r"\btor", tree_name="magnetics"
        )
        # Toroidal power supply takes time to turn on, from ~ -1.8 and should be
        # on by t=-1. So pick the time before that to calculate baseline
        baseline_indices = np.where(t_mag <= -1.8)
        btor = btor - np.mean(btor[baseline_indices])
        path = r"\mag_bp_coils.signals."
        # For each sensor:
        # 1. Subtract baseline offset
        # 2. Subtract btor pickup
        # 3. Interpolate bp onto shot timebase

        for i, bp13_name in enumerate(bp13_names):
            signal = params.mds_conn.get_data(path + bp13_name, tree_name="magnetics")
            if len(signal) == 1:
                raise CalculationError(f"No data for {bp13_name}")

            baseline = np.mean(signal[baseline_indices])
            signal = signal - baseline
            signal = signal - bp13_btor_pickup_coeffs[i] * btor
            bp13_signals[:, i] = interp1(t_mag, signal, params.times)

        # TODO: Examine edge case behavior of sign
        polarity = np.sign(np.mean(btor))
        btor_magnitude = btor * polarity
        btor_magnitude = interp1(t_mag, btor_magnitude, params.times)
        btor = interp1(t_mag, btor, params.times)  # Interpolate BT with sign

        # Create the 'design' matrix ('A') for the linear system of equations:
        # Bp(phi) = A1 + A2*sin(phi) + A3*cos(phi)
        ncoeffs = 3
        a = np.empty((len(bp13_names), ncoeffs))
        a[:, 0] = np.ones(4)
        a[:, 1] = np.sin(bp13_phi * np.pi / 180.0)
        a[:, 2] = np.cos(bp13_phi * np.pi / 180.0)
        coeffs = np.linalg.pinv(a) @ bp13_signals.T
        # The n=1 amplitude at each time is sqrt(A2^2 + A3^2)
        # The n=1 phase at each time is arctan(-A2/A3), using complex number
        # phasor formalism, exp(i(phi - delta))
        n_equal_1_amplitude = np.sqrt(coeffs[1, :] ** 2 + coeffs[2, :] ** 2)
        # TODO: Confirm arctan2 = atan2
        n_equal_1_phase = np.arctan2(-coeffs[1, :], coeffs[2, :])
        n_equal_1_normalized = n_equal_1_amplitude / btor_magnitude
        # INFO: Debugging purpose block of code at end of matlab file
        # INFO: n_equal_1_amplitude vs n_equal_1_mode
        output = {
            "n_equal_1_mode": n_equal_1_amplitude,
            "n_equal_1_normalized": n_equal_1_normalized,
            "n_equal_1_phase": n_equal_1_phase,
            "bt": btor,
        }
        return output

    @staticmethod
    def _get_densities(times, n_e, t_n, ip, t_ip, a_minor, t_a):
        """
        Calculate electron density, its time derivative, and the Greenwald fraction.

        Parameters
        ----------
        times : array_like
            Time points at which to interpolate the densities.
        n_e : array_like
            Electron density values.
        t_n : array_like
            Corresponding time values for electron density.
        ip : array_like
            Plasma current values.
        t_ip : array_like
            Corresponding time values for plasma current.
        a_minor : array_like
            Minor radius values.
        t_a : array_like
            Corresponding time values for minor radius.

        Returns
        -------
        dict
            A dictionary containing interpolated electron density (`n_e`),
            its time derivative (`dn_dt`), and the Greenwald fraction (`greenwald_fraction`).
        """
        if len(n_e) != len(t_n):
            raise CalculationError("n_e and t_n are different lengths")
        # get the gradient of n_E
        dn_dt = np.gradient(n_e, t_n)
        n_e = interp1(t_n, n_e, times)
        dn_dt = interp1(t_n, dn_dt, times)
        ip = -ip / 1e6  # Convert from A to MA and take positive value
        ip = interp1(t_ip, ip, times)
        a_minor = interp1(t_a, a_minor, times, bounds_error=False, fill_value=np.nan)
        # make sure aminor is not 0 or less than 0
        a_minor[a_minor <= 0] = 0.001
        n_g = abs(ip) / (np.pi * a_minor**2) * 1e20  # Greenwald density in m ^-3
        g_f = n_e / n_g
        output = {"n_e": n_e, "dn_dt": dn_dt, "greenwald_fraction": g_f}
        return output

    @staticmethod
    @physics_method(
        columns=["n_e", "dn_dt", "greenwald_fraction"],
        tokamak=Tokamak.CMOD,
    )
    def get_densities(params: PhysicsMethodParams):
        """
        Retrieve and calculate electron density and related parameters.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing electron density (`n_e`), its gradient (`dn_dt`),
            and the Greenwald fraction (`greenwald_fraction`).
        """
        # Line-integrated density
        n_e, t_n = params.mds_conn.get_data_with_dims(
            r".tci.results:nl_04", tree_name="electrons", astype="float64"
        )
        # Divide by chord length of ~0.6m to get line averaged density.
        # For future refernce, chord length is stored in
        # .01*\analysis::efit_aeqdsk:rco2v[3,*]
        n_e = np.squeeze(n_e) / 0.6
        ip, t_ip = params.mds_conn.get_data_with_dims(
            r"\ip", tree_name="magnetics", astype="float64"
        )
        a_minor, t_a = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:aminor", tree_name="_efit_tree", astype="float64"
        )

        output = CmodPhysicsMethods._get_densities(
            params.times, n_e, t_n, ip, t_ip, a_minor, t_a
        )
        return output

    @staticmethod
    def _get_efc_current(times, iefc, t_iefc):
        """
        Interpolate EFC current values at specified times.

        Parameters
        ----------
        times : array_like
            Time points at which to interpolate the EFC current.
        iefc : array_like
            EFC current values.
        t_iefc : array_like
            Corresponding time values for EFC current.

        Returns
        -------
        dict
            A dictionary containing interpolated EFC current.
        """
        output = {"i_efc": interp1(t_iefc, iefc, times, "linear")}
        return output

    @staticmethod
    @physics_method(columns=["i_efc"], tokamak=Tokamak.CMOD)
    def get_efc_current(params: PhysicsMethodParams):
        """
        Retrieve the error field correction (EFC) current for a given shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the EFC current (`i_efc`).
        """
        iefc, t_iefc = params.mds_conn.get_data_with_dims(
            r"\efc:u_bus_r_cur", tree_name="engineering"
        )
        output = CmodPhysicsMethods._get_efc_current(params.times, iefc, t_iefc)
        return output

    @staticmethod
    def _get_ts_parameters(times, ts_data, ts_time, ts_z, z_sorted=False):
        """
        Calculate the Thomson scattering temperature width parameters.

        Parameters
        ----------
        times : array_like
            Time points at which to interpolate the temperature width.
        ts_data : array_like
            2D array of Thomson scattering temperature data.
        ts_time : array_like
            Corresponding time values for the temperature data.
        ts_z : array_like
            Vertical coordinate values corresponding to the temperature data.
        z_sorted : bool, optional
            If True, assumes `ts_z` is already sorted. Default is False.

        Returns
        -------
        dict
            A dictionary containing the temperature width (`te_width`).
        """
        # sort z array
        if not z_sorted:
            idx = np.argsort(ts_z)
            ts_z = ts_z[idx]
            ts_data = ts_data[idx]
        # init output
        te_hwm = np.full(len(ts_time), np.nan)
        # select valid times
        (valid_times,) = np.where(ts_time > 0)
        # zero out nan values
        ts_data = np.nan_to_num(ts_data, copy=False, nan=0)
        # for each valid time
        for idx in valid_times:
            # select non-zero indices
            y = ts_data[:, idx]
            (ok_indices,) = np.where(y != 0)
            # skip if not enough points
            if len(ok_indices) < 3:
                continue
            # working arrays
            y = y[ok_indices]
            z = ts_z[ok_indices]
            # initial guess
            i = y.argmax()
            guess = [y[i], z[i], (z.max() - z.min()) / 3]
            # actual fit
            try:
                _, _, psigma = gaussian_fit(z, y, guess)
            except RuntimeError as exc:
                if str(exc).startswith("Optimal parameters not found"):
                    continue
                raise exc
            # store output
            te_hwm[idx] = np.abs(psigma)
        # rescale from sigma to HWHM
        # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        te_hwm *= np.sqrt(2 * np.log(2))
        # time interpolation
        te_hwm = interp1(ts_time, te_hwm, times)
        return {"te_width": te_hwm}

    @staticmethod
    @physics_method(columns=["te_width"], tokamak=Tokamak.CMOD)
    def get_ts_parameters(params: PhysicsMethodParams):
        """
        Retrieve Thomson scattering temperature width parameters.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the temperature width (`te_width`).
        """
        # TODO: Gaussian vs parabolic fit for te profile

        # Read in Thomson core temperature data, which is a 2-D array, with the
        # dependent dimensions being time and z (vertical coordinate)
        node_path = ".yag_new.results.profiles"

        ts_data, ts_time = params.mds_conn.get_data_with_dims(
            node_path + ":te_rz", tree_name="electrons"
        )
        ts_z = params.mds_conn.get_data(node_path + ":z_sorted", tree_name="electrons")

        output = CmodPhysicsMethods._get_ts_parameters(
            params.times, ts_data, ts_time, ts_z
        )
        return output

    @staticmethod
    def _get_peaking_factors(times, ts_time, ts_te, ts_ne, ts_z, efit_time, bminor, z0):
        """
        Calculate Te, ne, and pressure peaking factors given Thomson Scattering
        Te and ne measurements.

        Because the TS chords have uneven spacings, measurements are first interpolated
        to an array of equally spaced vertical positions and then used to calculate
        the peaking factors.

        Parameters:
        ----------
        times : array_like
            Requested time basis
        ts_time : array_like
            Time basis of the Thomson Scattering diagnostic
        ts_te : array_like
            Core and edge Te measurements from TS
        ts_ne : array_like
            Core and edge ne measurements from TS
        ts_z : array_like
            Vertical position of the core and edge TS chords
        efit_time : array_like
            Time basis of '_efit_tree'
        bminor : array_like
            Vertical minor radius from EFIT
        z0 : array_like
            Vertical position of the magnetic axis from EFIT

        Returns:
        ----------
        DataFrame of ne_peaking, Te_peaking, and pressure_peaking

        References:
        ----------
        - https://github.com/MIT-PSFC/disruption-py/blob/matlab/CMOD/matlab-core/get_peaking_factor_cmod.m  # pylint: disable=line-too-long
        - https://github.com/MIT-PSFC/disruption-py/issues/210
        - https://github.com/MIT-PSFC/disruption-py/pull/216
        - https://github.com/MIT-PSFC/disruption-py/pull/268

        Last major update by: William Wei on 8/19/2024

        """
        # Calculate ts_pressure
        ts_pressure = ts_te * ts_ne * 1.38e-23
        # Interpolate EFIT signals to TS time basis
        bminor = interp1(efit_time, bminor, ts_time)
        z0 = interp1(efit_time, z0, ts_time)

        # Calculate Te, ne, & pressure peaking factors
        te_pf = np.full(len(ts_time), np.nan)
        ne_pf = np.full(len(ts_time), np.nan)
        pressure_pf = np.full(len(ts_time), np.nan)
        (itimes,) = np.where((ts_time > 0) & (ts_time < times[-1]))
        for itime in itimes:
            ts_te_arr = ts_te[:, itime]
            ts_ne_arr = ts_ne[:, itime]
            ts_pressure_arr = ts_pressure[:, itime]
            # This gives identical results using either ts_te_arr or ts_ne_arr
            (indx,) = np.where(ts_ne_arr > 0)
            if len(indx) < 10:
                continue
            ts_te_arr = ts_te_arr[indx]
            ts_ne_arr = ts_ne_arr[indx]
            ts_pressure_arr = ts_pressure_arr[indx]
            ts_z_arr = ts_z[indx]
            sorted_indx = np.argsort(ts_z_arr)
            ts_z_arr = ts_z_arr[sorted_indx]
            ts_te_arr = ts_te_arr[sorted_indx]
            ts_ne_arr = ts_ne_arr[sorted_indx]
            ts_pressure_arr = ts_pressure_arr[sorted_indx]
            # Create equal-spacing array of ts_z_arr and interpolate TS profile on it
            # Skip if there's no EFIT zmagx data
            if np.isnan(z0[itime]):
                continue
            z_arr_equal_spacing = np.linspace(z0[itime], ts_z_arr[-1], len(ts_z_arr))
            te_arr_equal_spacing = interp1(ts_z_arr, ts_te_arr, z_arr_equal_spacing)
            ne_arr_equal_spacing = interp1(ts_z_arr, ts_ne_arr, z_arr_equal_spacing)
            pressure_arr_equal_spacing = interp1(
                ts_z_arr, ts_pressure_arr, z_arr_equal_spacing
            )
            # Calculate peaking factors
            (core_index,) = np.where(
                np.array((z_arr_equal_spacing - z0[itime]) < 0.2 * abs(bminor[itime]))
            )
            if len(core_index) < 2:
                continue
            te_pf[itime] = np.mean(te_arr_equal_spacing[core_index]) / np.mean(
                te_arr_equal_spacing
            )
            ne_pf[itime] = np.mean(ne_arr_equal_spacing[core_index]) / np.mean(
                ne_arr_equal_spacing
            )
            pressure_pf[itime] = np.mean(
                pressure_arr_equal_spacing[core_index]
            ) / np.mean(pressure_arr_equal_spacing)

        # Interpolate peaking factors to the requested time basis
        ne_pf = interp1(ts_time, ne_pf, times, "linear")
        te_pf = interp1(ts_time, te_pf, times, "linear")
        pressure_pf = interp1(ts_time, pressure_pf, times, "linear")
        return {
            "ne_peaking": ne_pf,
            "te_peaking": te_pf,
            "pressure_peaking": pressure_pf,
        }

    @staticmethod
    @physics_method(
        columns=["ne_peaking", "te_peaking", "pressure_peaking"],
        tokamak=Tokamak.CMOD,
    )
    def get_peaking_factors(params: PhysicsMethodParams):
        """
        Calculate peaking factors for electron density, electron temperature, and
        pressure.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing peaking factors for electron density (`ne_peaking`),
            temperature (`te_peaking`), and pressure (`pressure_peaking`).
        """
        use_ts_tci_calibration = False
        # Ignore shots on the blacklist
        if CmodPhysicsMethods.is_on_blacklist(params.shot_id):
            raise CalculationError("Shot is on blacklist")
        # Fetch data
        # Get EFIT geometry data
        z0 = 0.01 * params.mds_conn.get_data(
            r"\efit_aeqdsk:zmagx", tree_name="_efit_tree"
        )
        kappa = params.mds_conn.get_data(r"\efit_aeqdsk:kappa", tree_name="_efit_tree")
        aminor, efit_time = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
        )
        bminor = aminor * kappa

        # Get Te data and TS time basis
        node_ext = ".yag_new.results.profiles"
        ts_te_core, ts_time = params.mds_conn.get_data_with_dims(
            f"{node_ext}:te_rz", tree_name="electrons"
        )
        ts_te_core = ts_te_core * 1000  # [keV] -> [eV]
        ts_te_edge = params.mds_conn.get_data(r"\ts_te")
        ts_te = np.concatenate((ts_te_core, ts_te_edge)) * 11600  # [eV] -> [K]

        # Get ne data
        ts_ne_core = params.mds_conn.get_data(
            f"{node_ext}:ne_rz", tree_name="electrons"
        )
        ts_ne_edge = params.mds_conn.get_data(r"\ts_ne")
        ts_ne = np.concatenate((ts_ne_core, ts_ne_edge))

        # Get TS chord positions
        ts_z_core = params.mds_conn.get_data(
            f"{node_ext}:z_sorted", tree_name="electrons"
        )
        ts_z_edge = params.mds_conn.get_data(r"\fiber_z", tree_name="electrons")
        ts_z = np.concatenate((ts_z_core, ts_z_edge))
        # Make sure that there are equal numbers of edge position and edge temperature points
        if len(ts_z_edge) != ts_te_edge.shape[0]:
            raise CalculationError(
                "TS edge data and z positions are not the same length for shot"
            )

        # Calibrate ts_ne using TCI -- slow
        if use_ts_tci_calibration:
            # This shouldn't affect ne_PF (except if calib is not between 0.5 & 1.5)
            # because we're just multiplying ne by a constant
            (nl_ts1, nl_ts2, nl_tci1, nl_tci2, _, _) = (
                CmodThomsonDensityMeasure.compare_ts_tci(params)
            )
            if np.mean(nl_ts1) != 1e32 and np.mean(nl_ts2) != 1e32:
                nl_tci = np.concatenate((nl_tci1, nl_tci2))
                nl_ts = np.concatenate((nl_ts1 + nl_ts2))
                calib = np.mean(nl_tci) / np.mean(nl_ts)
            elif np.mean(nl_ts1) != 1e32 and np.mean(nl_ts2) == 1e32:
                calib = np.mean(nl_tci1) / np.mean(nl_ts1)
            elif np.mean(nl_ts1) == 1e32 and np.mean(nl_ts2) != 1e32:
                calib = np.mean(nl_tci2) / np.mean(nl_ts2)
            else:
                calib = np.nan

            if 0.5 < calib < 1.5:
                ts_ne *= calib
            else:
                raise CalculationError(
                    "Density calibration error exceeds acceptable range"
                )

        return CmodPhysicsMethods._get_peaking_factors(
            params.times, ts_time, ts_te, ts_ne, ts_z, efit_time, bminor, z0
        )

    @staticmethod
    def _get_te_profile_params_ece(
        times,
        gpc1_te_data,
        gpc1_te_time,
        gpc1_rad_data,
        gpc1_rad_time,
        gpc2_te_data,
        gpc2_te_time,
        gpc2_rad_data,
        gpc2_rad_time,
        efit_time,
        r0,
        aminor,
        btor,
        t_mag,
        lh_power,
        lh_time,
    ):
        """
        Calculates Te PF and width from ECE data using the two GPC diagnostic systems.
        GPC diagnostics look at the mid-plane, and each channel detects a different
        emitted frequency associated with the second harmonic, which depends on B and
        therefore R.
        - te_width is the half-width at half-max of a Gaussian fit of the Te profile
        - te_core_vs_avg is defined as mean(core)/mean(all) where core bins are defined
          as those w/ |R - R0| < 0.2*a of the magnetic axis.
        - te_edge_vs_avg is defined as mean(edge)/mean(all) where edge bins are defined as
          those with 0.8*a < |R - R0| < a
        For core and edge vs. average calculations, different shots can have different
        radial sampling, and during a few experiments on C-Mod, Bt was changed during
        the shot, changing the radial sampling. Different radial samplings can have
        different proportions of core to edge sampling, which affects the mean Te over
        the whole profile, biasing the core vs average and edge vs average statistics.
        Therefore, we use a uniformly sampled radial basis from R0 to R0+a. We use many
        interpolated radial points to minimize artifacts caused by a point moving
        across the arbitrary core or edge boundary.

        ECE as a Te profile diagnostic can suffer from several artifacts:
        Artifacts currently NOT explicitly checked for
        - Density cutoffs: High ne plasmas (typically H-modes) can have an ECE cutoff.
          According to Amanda Hubbard, "what you wil see is a section of profile which
          is much LOWER than Thomson Scattering, for some portion of the LFS profile
          (typically starting around r/a 0.8?). In this case ECE cannot be used." An
          example shot with ECE cutoffs is 1140226024 (Calibration of Thomson density
          using ECE cutoffs). Because the critical density is proportional to B^2,
          shots with B = 5.4 T on axis would need to have very high densities to
          experience a cutoff in the profile. We could look for cutoffs by comparing
          the B profile to the ne profile and checking that ne < ncrit throughout the
          profile; however, a simpler check for now is to ignore shots with B < 4.5 T
          and assume there are no cutoffs with B >= 4.5 T.
        Artifacts currently checked for
        - Non-aligned grating: The gratings were usually aligned for radial coverage
          assuming Bt=5.4T. For low Bt shots (like 2.8T), sometimes the gratings were
          adjusted, sometimes not. Low Bt shots also tend to have low signal and often
          experience density cutoffs. Therefore, ECE should be avoided in automated
          calculations for low Bt shots.
        - Non-thermal emission. The calculation of Te vs. r assumes that the second
          harmonic emission can be modeled as black-body emission, which assumes the
          electrons are in thermal equilibrium. On C-Mod, non-thermal emission results
          in an apparent Te that goes UP towards the edge and in the SOL, which is
          actually downshifted non-thermal emission from deeper in the core.
          Significant runaway populations and LHCD lead to non-thermal artifacts.
          Occasionally low ne shots also had non-thermal artifacts.
        - Harmonic overlap: Certain channels can pick up emission from different
          harmonics from other regions of the plasma. Generally channels with R < 0.6 m
          suffer from overlap with 3rd harmonic emission from the core. This leads to
          an apparently higher Te for R < 0.6 m than in reality. The gratings were
          usually aligned to measure the profile from the core outwards for this
          reason.

        Parameters
        ----------
        times : array_like
            Requested time basis
        gpc1_te_array: array_like
            Te measurements from GPC diagnostic
        gpc1_te_time: array_like
            Time basis of GPC Te measurements
        gpc1_rad_data: array_like
            Radial positions corresponding to GPC channels
        gpc1_rad_time: array_like
            Time basis of GPC channel radial positions
        gpc2_te_array: array_like
            Te measurements from GPC2 diagnostic
        gpc2_te_time: array_like
            Time basis of GPC2 Te measurements
        gpc2_rad_data: array_like
            Radial positions corresponding to GPC2 channels
        gpc2_rad_time: array_like
            Time basis of GPC2 channel radial positions
        efit_time : array_like
            Time basis of '_efit_tree'
        r0 : array_like
            Radial position of the magnetic axis from EFIT
        aminor : array_like
            Horizontal minor radius from EFIT
        btor: array_like
            On-axis toroidal field from magnetics
        t_mag: array_like
            Time basis of magnetic diagnostic measurements
        lh_power: array_like
            Lower hybrid power
        lh_time: array_like
            Time basis of lower hybrid power

        Returns
        -------
        Dictionary of ne_peaking, Te_peaking, and pressure_peaking

        Sources:
        - https://github.com/MIT-PSFC/disruption-py/blob/matlab/CMOD/matlab-core/
          get_ECE_data_cmod.m
        - K. Zhurovich, et. al. "Calibration of Thomson scattering systems using
          electron cyclotron emission cutoff data," Rev. Sci. Instrum., vol. 76, no. 5,
          p. 053506, 2005, doi: 10.1063/1.1899311.
        - https://github.com/MIT-PSFC/disruption-py/pull/260

        Last Major Update: Henry Wietfeldt (08/28/24), (PR: #260)
        """

        # Constants
        core_bound_factor = 0.2
        edge_bound_factor = 0.8
        min_okay_channels = 9
        min_te = 0.02  # [keV]
        min_btor = 4.5  # [T]
        max_lh_power = 1.0  # [kW]
        min_r_to_avoid_harmonic_overlap = 0.6  # [m]
        rising_tail_factor = 1.2

        # Only use EFITs starting after the GPC diagnostic has profiles.
        if len(gpc1_rad_time) > 0:
            efit_time = efit_time[
                efit_time >= max(np.max(gpc1_rad_time[:, 0]), gpc2_rad_time[0])
            ]
        else:
            efit_time = efit_time[efit_time >= gpc2_rad_time[0]]

        # Interpolate GPC data onto efit timebase. Timebase for radial measurements is
        # slower than efit but radial positions are approx. stable so linear
        # interpolation is safe.
        n_channels = gpc1_te_data.shape[0]
        gpc1_te = np.full((n_channels, len(efit_time)), np.nan)
        gpc1_rad = np.full((n_channels, len(efit_time)), np.nan)
        for i in range(n_channels):
            gpc1_te[i, :] = interp1(gpc1_te_time[i, :], gpc1_te_data[i, :], efit_time)
            if len(gpc1_rad_data[i, :]) > 1:
                gpc1_rad[i, :] = interp1(
                    gpc1_rad_time[i, :], gpc1_rad_data[i, :], efit_time
                )
            else:
                gpc1_rad[i, :] = np.full(len(efit_time), np.nan)

        n_channels = gpc2_te_data.shape[0]
        gpc2_te = np.full((n_channels, len(efit_time)), np.nan)
        gpc2_rad = np.full((n_channels, len(efit_time)), np.nan)
        for i in range(n_channels):
            gpc2_te[i, :] = interp1(gpc2_te_time, gpc2_te_data[i, :], efit_time)
            if len(gpc2_rad_data[i, :]) > 1:
                gpc2_rad[i, :] = interp1(gpc2_rad_time, gpc2_rad_data[i, :], efit_time)
            else:
                gpc2_rad[i, :] = np.full(len(efit_time), np.nan)

        # Combine GPC systems and extend the last radii measurement up until the last
        # EFIT. Radii depend on Bt, which should be stable until the current quench.
        te = np.concatenate((gpc1_te, gpc2_te), axis=0)
        radii = np.concatenate((gpc1_rad, gpc2_rad), axis=0)
        indx_last_rad = np.argmax(efit_time > gpc2_rad_time[-1]) - 1
        for i in range(len(radii)):
            radii[i, indx_last_rad + 1 :] = radii[i, indx_last_rad]

        # Remaining calculations loop over time then radii so transpose for efficient
        # caching
        te = te.T
        radii = radii.T
        for i in range(len(efit_time)):
            sorted_index = np.argsort(radii[i, :])
            radii[i, :] = radii[i, sorted_index]
            te[i, :] = te[i, sorted_index]

        # Time slices with low Btor are unreliable because gratings are often not
        # aligned to field, signal is low, and there are frequent density cutoffs.
        # Time slices with LH heating are unreliable because direct electron heating
        # leads to non-thermal emission
        btor = interp1(t_mag, btor, efit_time)
        if len(lh_time) > 1:
            lh_power = interp1(lh_time, lh_power, efit_time)
        else:
            lh_power = np.zeros(len(efit_time))
        lh_power = np.nan_to_num(lh_power, nan=0.0)
        (okay_time_indices,) = np.where(
            (np.abs(btor) > min_btor) & (lh_power < max_lh_power)
        )

        # Main loop for calculations
        te_core_vs_avg = np.full(len(efit_time), np.nan)
        te_edge_vs_avg = np.full(len(efit_time), np.nan)
        te_hwhm = np.full(len(efit_time), np.nan)
        for i in okay_time_indices:
            # Only consider points that are likely to accurately measure Te
            calib_indices = (te[i, :] > min_te) & (radii[i, :] > 0)
            harmonic_overlap_indices = radii[i, :] < min_r_to_avoid_harmonic_overlap
            nonthermal_overlap_indices = np.full(len(radii[i, :]), False)

            # Identify rising tail (overlap with non-thermal emission). Finding the min
            # Te near the edge and checking outwards for a rising tail seems to do well
            calib_edge = calib_indices & (
                radii[i, :] > r0[i] + edge_bound_factor * aminor[i]
            )
            if np.sum(calib_edge) > 0:
                te_edge = np.min(te[i, calib_edge])
                indx_edge = np.argmin(np.abs(te[i, :] - te_edge))
                for j in range(len(te[i, :]) - 1 - indx_edge):
                    if te[i, indx_edge + j + 1] > rising_tail_factor * te[i, indx_edge]:
                        nonthermal_overlap_indices[indx_edge + j + 1] = True
            okay_indices = (
                calib_indices
                & (~harmonic_overlap_indices)
                & (~nonthermal_overlap_indices)
            )

            if np.sum(okay_indices) > min_okay_channels:
                # Estimate Te width using Gaussian fit with center fixed on mag. axis
                r = radii[i, okay_indices]
                y = te[i, okay_indices]
                guess = [y.max(), (y.max() - y.min()) / 3]
                try:
                    pmu = r0[i]
                    _, psigma = gaussian_fit_with_fixed_mean(pmu, r, y, guess)
                except RuntimeError as exc:
                    if str(exc).startswith("Optimal parameters not found"):
                        continue
                    raise exc

                # rescale from sigma to HWHM
                # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
                te_hwhm[i] = np.abs(psigma) * np.sqrt(2 * np.log(2))

                # Calculate core/edge vs. average using uniformly sampled radial basis
                r_equal_spaced = np.linspace(r0[i], r0[i] + aminor[i], 100)
                te_equal_spaced = interp1(
                    r, y, r_equal_spaced, fill_value=(y[0], y[-1])
                )
                core_indices = (
                    np.abs(r_equal_spaced - r0[i]) < core_bound_factor * aminor[i]
                ) & (~np.isnan(te_equal_spaced))
                edge_indices = (
                    np.abs(r_equal_spaced - r0[i]) > edge_bound_factor * aminor[i]
                ) & (~np.isnan(te_equal_spaced))
                if np.sum(core_indices) > 0:
                    te_core_vs_avg[i] = np.nanmean(
                        te_equal_spaced[core_indices]
                    ) / np.nanmean(te_equal_spaced)
                if np.sum(edge_indices) > 0:
                    te_edge_vs_avg[i] = np.nanmean(
                        te_equal_spaced[edge_indices]
                    ) / np.nanmean(te_equal_spaced)

        te_core_vs_avg = interp1(efit_time, te_core_vs_avg, times)
        te_edge_vs_avg = interp1(efit_time, te_edge_vs_avg, times)
        te_hwhm = interp1(efit_time, te_hwhm, times)
        return {
            "te_core_vs_avg_ece": te_core_vs_avg,
            "te_edge_vs_avg_ece": te_edge_vs_avg,
            "te_width_ece": te_hwhm,
        }

    @staticmethod
    @physics_method(
        columns=["te_core_vs_avg_ece", "te_edge_vs_avg_ece", "te_width_ece"],
        tokamak=Tokamak.CMOD,
    )
    def get_te_profile_params_ece(params: PhysicsMethodParams):
        """
        Gets MDSplus data to be used in the calculations of te profile parameters
        from ECE data
        Parameters
        ----------
        params: PhysicsMethodParams
            The parameters storing the requested time base, shot id, etc
        Returns
        ----------
        Output of get_te_profile_params_ece(), which processes the MDSplus data

        Last Major Update: Henry Wietfeldt (8/28/24)
        """

        # Constants
        n_gpc1_channels = 9

        # Get magnetic axis data from EFIT
        r0 = 0.01 * params.mds_conn.get_data(
            r"\efit_aeqdsk:rmagx", tree_name="_efit_tree"
        )  # [cm] -> [m]
        aminor, efit_time = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
        )  # [m], [s]

        # Btor and LH Power used for filtering okay time slices
        btor, t_mag = params.mds_conn.get_data_with_dims(
            r"\btor", tree_name="magnetics"
        )
        try:
            lh_power, lh_time = params.mds_conn.get_data_with_dims(
                ".results:netpow", tree_name="lh"
            )  # [kW], [s]
        except mdsExceptions.MdsException:
            # When LH power is off, it's often not stored in tree or it's a single 0.
            lh_power = 0.0
        if not isinstance(lh_power, np.ndarray):
            lh_time = np.copy(efit_time)
            lh_power = np.zeros(len(efit_time))

        # Read in Te profile measurements from 9 GPC1 ("GPC" in MDSplus tree) channels
        node_path = ".ece.gpc_results"
        gpc1_te_data = []
        gpc1_te_time = []
        gpc1_rad_data = []
        gpc1_rad_time = []
        for i in range(n_gpc1_channels):
            try:
                te_data, te_time = params.mds_conn.get_data_with_dims(
                    node_path + ".te:te" + str(i + 1), tree_name="electrons"
                )  # [keV], [s]
                rad_data, rad_time = params.mds_conn.get_data_with_dims(
                    node_path + ".rad:r" + str(i + 1), tree_name="electrons"
                )  # [m], [s]
                # For C-Mod shot 1120522025 (and maybe others), rad_time is strings.
                # Don't use channel in that case
                if np.issubdtype(rad_time.dtype, np.floating):
                    gpc1_te_data.append(te_data)
                    gpc1_te_time.append(te_time)
                    gpc1_rad_data.append(rad_data)
                    gpc1_rad_time.append(rad_time)
            except mdsExceptions.MdsException:
                continue
        gpc1_te_data = np.array(gpc1_te_data)
        gpc1_te_time = np.array(gpc1_te_time)
        gpc1_rad_data = np.array(gpc1_rad_data)
        gpc1_rad_time = np.array(gpc1_rad_time)

        # Read in Te profile measurements from GPC2 (19 channels)
        node_path = ".gpc_2.results"
        gpc2_te_data, gpc2_te_time = params.mds_conn.get_data_with_dims(
            node_path + ":gpc2_te", tree_name="electrons"
        )  # [keV], [s]
        gpc2_rad_data, gpc2_rad_time = params.mds_conn.get_data_with_dims(
            node_path + ":radii", tree_name="electrons"
        )  # [m], [s]

        return CmodPhysicsMethods._get_te_profile_params_ece(
            params.times,
            gpc1_te_data,
            gpc1_te_time,
            gpc1_rad_data,
            gpc1_rad_time,
            gpc2_te_data,
            gpc2_te_time,
            gpc2_rad_data,
            gpc2_rad_time,
            efit_time,
            r0,
            aminor,
            btor,
            t_mag,
            lh_power,
            lh_time,
        )

    @staticmethod
    @physics_method(
        columns=["prad_peaking"],
        tokamak=Tokamak.CMOD,
    )
    def get_prad_peaking(params: PhysicsMethodParams):
        """
        Calculate the peaking factor for radiated power.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the peaking factor for radiated power (`prad_peaking`).
        """
        prad_peaking = np.full(len(params.times), np.nan)
        nan_output = {"prad_peaking": prad_peaking}
        r0 = 0.01 * params.mds_conn.get_data(
            r"\efit_aeqdsk:rmagx", tree_name="_efit_tree"
        )
        z0 = 0.01 * params.mds_conn.get_data(
            r"\efit_aeqdsk:zmagx", tree_name="_efit_tree"
        )
        aminor, efit_time = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
        )
        got_axa = False
        try:
            bright_axa, t_axa, r_axa = params.mds_conn.get_data_with_dims(
                r"\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT",
                tree_name="spectroscopy",
                dim_nums=[1, 0],
            )
            z_axa = params.mds_conn.get_data(
                r"\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:Z_O",
                tree_name="spectroscopy",
            )
            good_axa = params.mds_conn.get_data(
                r"\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:GOOD",
                tree_name="spectroscopy",
            )
            got_axa = True
        except mdsExceptions.MdsException:
            params.logger.debug("[Shot %s]: Failed to get AXA data", params.shot_id)
        got_axj = False
        try:
            bright_axj, t_axj, r_axj = params.mds_conn.get_data_with_dims(
                r"\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXJ:BRIGHT",
                tree_name="spectroscopy",
                dim_nums=[1, 0],
            )
            z_axj = params.mds_conn.get_data(
                r"\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXJ:Z_O",
                tree_name="spectroscopy",
            )
            good_axj = params.mds_conn.get_data(
                r"\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXJ:GOOD",
                tree_name="spectroscopy",
            )
            got_axj = True
        except mdsExceptions.MdsException:
            params.logger.debug("[Shot %s]: Failed to get AXJ data", params.shot_id)
        if not (got_axa or got_axj):
            return nan_output
        a_minor = interp1(efit_time, aminor, params.times)
        r0 = interp1(efit_time, r0, params.times)
        z0 = interp1(efit_time, z0, params.times)
        if got_axa:
            good_axa = np.where(good_axa > 0)[0]
            bright_axa = bright_axa[:, good_axa]
            axa_interp = np.full((bright_axa.shape[1], len(params.times)), np.nan)
            r_axa = r_axa[good_axa]
            for i in range(bright_axa.shape[1]):
                interped = interp1(t_axa.T, bright_axa[:, i], params.times.T)
                indx = np.where(interped < 0)
                interped[indx] = np.nan
                axa_interp[i, :] = interped
        if got_axj:
            good_axj = np.where(good_axj > 0)[0]
            bright_axj = bright_axj[:, good_axj]
            axj_interp = np.full((bright_axj.shape[1], len(params.times)), np.nan)
            r_axj = r_axj[good_axj]
            for i in range(bright_axj.shape[1]):
                interped = interp1(t_axj.T, bright_axj[:, i], params.times.T)
                indx = np.where(interped < 0)
                interped[indx] = np.nan
                axj_interp[i, :] = interped
        for i in range(len(params.times)):
            core_radiation = np.array([])
            all_radiation = np.array([])
            if got_axa:
                axa_dist = np.sqrt((r_axa - r0[i]) ** 2 + (z0[i] - z_axa) ** 2)
                axa_core_index = axa_dist < 0.2 * a_minor[i]
                core_radiation = np.append(
                    core_radiation, axa_interp[axa_core_index, i]
                )
                all_radiation = np.append(all_radiation, axa_interp[:, i])
            if got_axj:
                axj_dist = np.sqrt((r_axj - r0[i]) ** 2 + (z0[i] - z_axj) ** 2)
                axj_core_index = axj_dist < 0.2 * a_minor[i]
                core_radiation = np.append(
                    core_radiation, axj_interp[axj_core_index, i]
                )
                all_radiation = np.append(all_radiation, axj_interp[:, i])
            with warnings.catch_warnings():
                warnings.filterwarnings(action="ignore", message="Mean of empty slice")
                prad_peaking[i] = np.nanmean(core_radiation) / np.nanmean(all_radiation)
        return {"prad_peaking": prad_peaking}

    # TODO: get more accurate description of soft x-ray data
    @staticmethod
    @physics_method(columns=["sxr"], tokamak=Tokamak.CMOD)
    def get_sxr_data(params: PhysicsMethodParams):
        """
        Retrieve soft X-ray (SXR) data from array 1 chord 16 for a given shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDSplus connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the soft X-ray data (`sxr`).
        """
        sxr, t_sxr = params.mds_conn.get_data_with_dims(
            r"\top.brightnesses.array_1:chord_16",
            tree_name="xtomo",
            astype="float64",
        )
        sxr = interp1(t_sxr, sxr, params.times)
        return {"sxr": sxr}

    @staticmethod
    def is_on_blacklist(shot_id: int) -> bool:
        """
        TODO why will these shots cause `_get_peaking_factors`,
        `_get_peaking_factors_no_tci`, and `_get_edge_parameters` to fail?
        """
        if (
            1120000000 < shot_id < 1120213000
            or 1140000000 < shot_id < 1140227000
            or 1150000000 < shot_id < 1150610000
            or 1160000000 < shot_id < 1160303000
        ):
            return True
        return False
