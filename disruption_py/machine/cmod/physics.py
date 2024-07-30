#!/usr/bin/env python3

import traceback
import warnings
from importlib import resources

import numpy as np
import pandas as pd
import scipy as sp
from MDSplus import mdsExceptions

import disruption_py.data
from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import gaussian_fit, interp1, smooth
from disruption_py.machine.tokamak import Tokamak

warnings.filterwarnings("error", category=RuntimeWarning)


class CmodPhysicsMethods:
    @staticmethod
    @cache_method
    def get_active_wire_segments(params: PhysicsMethodParams):
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
                    f'getnci($, "STATE")', arguments=node_path + ":SEG_NUM"
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
    def _get_time_until_disrupt(params: PhysicsMethodParams):
        time_until_disrupt = [np.nan]
        if params.disrupted:
            time_until_disrupt = params.disruption_time - params.times
        return {"time_until_disrupt": time_until_disrupt}

    @staticmethod
    def get_ip_parameters(times, ip, magtime, ip_prog, pcstime):
        """Calculates actual and programmed current as well as their derivatives
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
    def _get_ip_parameters(params: PhysicsMethodParams):
        # Automatically generated
        active_segments = CmodPhysicsMethods.get_active_wire_segments(params=params)

        # Default PCS timebase is 1 KHZ
        pcstime = np.array(np.arange(-4, 12.383, 0.001))
        ip_prog = np.full(pcstime.shape, np.nan)

        # For each activate segment:
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
                    except mdsExceptions.MdsException as e:
                        params.logger.warning(
                            [
                                f"[Shot {params.shot_id}]: Error getting PID gains"
                                + f" for wire {wire_index}"
                            ]
                        )
                        params.logger.debug(
                            [f"[Shot {params.shot_id}]: {traceback.format_exc()}"]
                        )
                    break  # Break out of wire_index loop
        ip, magtime = params.mds_conn.get_data_with_dims(
            r"\ip", tree_name="magnetics", astype="float64"
        )
        output = CmodPhysicsMethods.get_ip_parameters(
            params.times, ip, magtime, ip_prog, pcstime
        )
        return output

    @staticmethod
    def get_z_parameters(times, z_prog, pcstime, z_error_without_ip, ip, dpcstime):
        """Get values of Z_error, Z_prog, and derived signals from plasma control
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
    def _get_z_parameters(params: PhysicsMethodParams):
        pcstime = np.array(np.arange(-4, 12.383, 0.001))
        z_prog = np.empty(pcstime.shape)
        z_prog.fill(np.nan)
        z_prog_temp = z_prog.copy()
        z_wire_index = -1
        active_wire_segments = CmodPhysicsMethods.get_active_wire_segments(
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
                    except mdsExceptions.MdsException as e:
                        params.logger.debug(
                            f"[Shot {params.shot_id}]: {traceback.format_exc()}"
                        )
                        continue  # TODO: Consider raising appropriate error
                else:
                    continue
                break
        if z_wire_index == -1:
            # TODO: Make appropriate error
            raise ValueError("No ZCUR wire was found")
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
        for i in range(len(active_wire_segments)):
            segment, start = active_wire_segments[i]
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
        return CmodPhysicsMethods.get_z_parameters(
            params.times, z_prog, pcstime, z_error_without_ip, ip, dpcstime
        )

    @staticmethod
    def get_ohmic_parameters(
        times, v_loop, v_loop_time, li, efittime, dip_smoothed, ip
    ):
        """Calculate the ohmic power from the loop voltage, inductive voltage, and
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
            The inductance of the loop.
        efittime : array_like
            The times at which the inductance was measured.
        dip_smoothed : array_like
            The smoothed plasma current.
        ip : array_like
            The plasma current.

        Returns
        -------
        p_ohm : array_like
            The ohmic power.
        v_loop : array_like
            The loop voltage.

        Original Authors
        ----------------


        """
        # For simplicity, we use R0 = 0.68 m, but we could use \efit_aeqdsk:rmagx
        R0 = 0.68
        inductance = 4.0 * np.pi * 1.0e-7 * R0 * li / 2.0
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
    def _get_ohmic_parameters(params: PhysicsMethodParams):
        v_loop, v_loop_time = params.mds_conn.get_data_with_dims(
            r"\top.mflux:v0", tree_name="analysis", astype="float64"
        )
        if len(v_loop_time) <= 1:
            return {
                "p_oh": np.zeros(len(params.times)),
                "v_loop": np.zeros(len(params.times)),
            }
        li, efittime = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:li", tree_name="_efit_tree", astype="float64"
        )
        ip_parameters = CmodPhysicsMethods._get_ip_parameters(params=params)
        output = CmodPhysicsMethods.get_ohmic_parameters(
            params.times,
            v_loop,
            v_loop_time,
            li,
            efittime,
            ip_parameters["dip_smoothed"],
            ip_parameters["ip"],
        )
        return output

    @staticmethod
    def get_power(times, p_lh, t_lh, p_icrf, t_icrf, p_rad, t_rad, p_ohm):
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
    def _get_power(params: PhysicsMethodParams):
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
            except (mdsExceptions.TreeFOPENR, mdsExceptions.TreeNNF) as e:
                continue
        p_oh = CmodPhysicsMethods._get_ohmic_parameters(params=params)["p_oh"]
        output = CmodPhysicsMethods.get_power(params.times, *values, p_oh)
        return output

    @staticmethod
    def get_kappa_area(times, aminor, area, a_times):
        output = {"kappa_area": interp1(a_times, area / (np.pi * aminor**2), times)}
        return output

    @staticmethod
    @physics_method(columns=["kappa_area"], tokamak=Tokamak.CMOD)
    def _get_kappa_area(params: PhysicsMethodParams):
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
        output = CmodPhysicsMethods.get_kappa_area(params.times, aminor, area, times)
        return output

    @staticmethod
    def get_rotation_velocity(times, intensity, time, vel, hirextime):
        """
        Uses spectroscopy graphs of ionized(to hydrogen and helium levels) Argon
        to calculate velocity. Because of the heat profile of the plasma, suitable
        measurements are only found near the center
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

    # TODO: Calculate v_mid
    @staticmethod
    @physics_method(columns=["v_0"], tokamak=Tokamak.CMOD)
    def _get_rotation_velocity(params: PhysicsMethodParams):
        nan_output = {"v_0": [np.nan]}
        with resources.path(disruption_py.data, "lock_mode_calib_shots.txt") as fio:
            calibrated = pd.read_csv(fio)
        # Check to see if shot was done on a day where there was a locked
        # mode HIREX calibration by cross checking with list of calibrated
        # runs. If not calibrated, return NaN outputs.
        if params.shot_id not in calibrated:
            return nan_output
        try:
            intensity, time = params.mds_conn.get_data_with_dims(
                ".hirex_sr.analysis.a:int", tree_name="spectroscopy", astype="float64"
            )
            vel, hirextime = params.mds_conn.get_data_with_dims(
                ".hirex_sr.analysis.a:vel", tree_name="spectroscopy", astype="float64"
            )
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_id}]: Failed to open necessary tress for "
                + f"rotational velocity calculations."
            )
            params.logger.debug(f"[Shot {params.shot_id}]: {traceback.format_exc()}")
            return nan_output
        output = CmodPhysicsMethods.get_rotation_velocity(
            params.times, intensity, time, vel, hirextime
        )
        return output

    # TODO: Split into static and instance method
    @staticmethod
    def get_n_equal_1_amplitude():
        pass

    # TODO: Try catch failure to get BP13 sensors
    @staticmethod
    @physics_method(
        columns=["n_equal_1_mode", "n_equal_1_normalized", "n_equal_1_phase", "bt"],
        tokamak=Tokamak.CMOD,
    )
    def _get_n_equal_1_amplitude(params: PhysicsMethodParams):
        """Calculate n=1 amplitude and phase.

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

        # Only calculate n=1 amplitude if all sensors have data
        valid_sensors = True
        for i in range(len(bp13_names)):
            try:
                signal = params.mds_conn.get_data(
                    path + bp13_names[i], tree_name="magnetics"
                )
                if len(signal) == 1:
                    params.logger.warning(
                        f"[Shot {params.shot_id}] Only one data point for "
                        + f"{bp13_names[i]} Returning nans."
                    )
                    return {
                        "n_equal_1_mode": [np.nan],
                        "n_equal_1_normalized": [np.nan],
                        "n_equal_1_phase": [np.nan],
                        "bt": [np.nan],
                    }
                baseline = np.mean(signal[baseline_indices])
                signal = signal - baseline
                signal = signal - bp13_btor_pickup_coeffs[i] * btor
                bp13_signals[:, i] = interp1(t_mag, signal, params.times)
            except mdsExceptions.TreeNODATA as e:
                params.logger.warning(
                    f"[Shot {params.shot_id}] No data for {bp13_names[i]}"
                )
                params.logger.debug(f"[Shot {params.shot_id}] {e}")
                valid_sensors = False
        # TODO: Examine edge case behavior of sign
        polarity = np.sign(np.mean(btor))
        btor_magnitude = btor * polarity
        btor_magnitude = interp1(t_mag, btor_magnitude, params.times)
        btor = interp1(t_mag, btor, params.times)  # Interpolate BT with sign
        if valid_sensors:
            # Create the 'design' matrix ('A') for the linear system of equations:
            # Bp(phi) = A1 + A2*sin(phi) + A3*cos(phi)
            ncoeffs = 3
            A = np.empty((len(bp13_names), ncoeffs))
            A[:, 0] = np.ones(4)
            A[:, 1] = np.sin(bp13_phi * np.pi / 180.0)
            A[:, 2] = np.cos(bp13_phi * np.pi / 180.0)
            coeffs = np.linalg.pinv(A) @ bp13_signals.T
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
    def get_densities(times, n_e, t_n, ip, t_ip, a_minor, t_a):
        if len(n_e) != len(t_n):
            return {
                "n_e": [np.nan],
                "dn_dt": [np.nan],
                "greenwald_fraction": [np.nan],
            }
        # get the gradient of n_E
        dn_dt = np.gradient(n_e, t_n)
        n_e = interp1(t_n, n_e, times)
        dn_dt = interp1(t_n, dn_dt, times)
        ip = -ip / 1e6  # Convert from A to MA and take positive value
        ip = interp1(t_ip, ip, times)
        a_minor = interp1(t_a, a_minor, times, bounds_error=False, fill_value=np.nan)
        # make sure aminor is not 0 or less than 0
        a_minor[a_minor <= 0] = 0.001
        n_G = abs(ip) / (np.pi * a_minor**2) * 1e20  # Greenwald density in m ^-3
        g_f = n_e / n_G
        output = {"n_e": n_e, "dn_dt": dn_dt, "greenwald_fraction": g_f}
        return output

    @staticmethod
    @physics_method(
        columns=["n_e", "dn_dt", "greenwald_fraction"],
        tokamak=Tokamak.CMOD,
    )
    def _get_densities(params: PhysicsMethodParams):
        try:
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
        except Exception as e:
            params.logger.debug(f"[Shot {params.shot_id}] {e}")
            params.logger.warning(f"[Shot {params.shot_id}] No density data")
            # TODO: Handle this case
            raise NotImplementedError(
                "Can't currently handle failure of grabbing density data"
            )
        output = CmodPhysicsMethods.get_densities(
            params.times, n_e, t_n, ip, t_ip, a_minor, t_a
        )
        return output

    @staticmethod
    def get_efc_current(times, iefc, t_iefc):
        output = {"i_efc": interp1(t_iefc, iefc, times, "linear")}
        return output

    @staticmethod
    @physics_method(columns=["i_efc"], tokamak=Tokamak.CMOD)
    def _get_efc_current(params: PhysicsMethodParams):
        try:
            iefc, t_iefc = params.mds_conn.get_data_with_dims(
                r"\efc:u_bus_r_cur", tree_name="engineering"
            )
        except Exception as e:
            params.logger.debug(f"[Shot {params.shot_id}] {traceback.format_exc()}")
            return {"i_efc": [np.nan]}
        output = CmodPhysicsMethods.get_efc_current(params.times, iefc, t_iefc)
        return output

    @staticmethod
    def get_Ts_parameters(times, ts_data, ts_time, ts_z, z_sorted=False):
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
    def _get_Ts_parameters(params: PhysicsMethodParams):
        # TODO: Guassian vs parabolic fit for te profile

        # Read in Thomson core temperature data, which is a 2-D array, with the
        # dependent dimensions being time and z (vertical coordinate)
        node_path = ".yag_new.results.profiles"
        try:
            ts_data, ts_time = params.mds_conn.get_data_with_dims(
                node_path + ":te_rz", tree_name="electrons"
            )
            ts_z = params.mds_conn.get_data(
                node_path + ":z_sorted", tree_name="electrons"
            )
        except mdsExceptions.MdsException as e:
            params.logger.debug(f"[Shot {params.shot_id}] {traceback.format_exc()}")
            return {"te_width": [np.nan]}
        output = CmodPhysicsMethods.get_Ts_parameters(
            params.times, ts_data, ts_time, ts_z
        )
        return output

    @staticmethod
    def get_peaking_factors(times, TS_time, TS_Te, TS_ne, TS_z, efit_time, bminor, z0):
        """
        Calculate Te, ne, and pressure peaking factors given Thomson Scattering Te and ne measurements.

        Because the TS chords have uneven spacings, measurements are first interpolated to an array of
        equally spaced vertical positions and then used to calculate the peaking factors.

        Currently, only the Te_peaking feature has been implemented.

        Parameters:
        ----------
        times : array_like
            Requested time basis
        TS_time : array_like
            Time basis of the Thomson Scattering diagnostic
        TS_Te : array_like
            Core and edge Te measurements from TS
        TS_ne : array_like
            Core and edge ne measurements from TS
        TS_z : array_like
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
        - https://github.com/MIT-PSFC/disruption-py/blob/matlab/CMOD/matlab-core/get_peaking_factor_cmod.m
        - https://github.com/MIT-PSFC/disruption-py/issues/210
        - https://github.com/MIT-PSFC/disruption-py/pull/216

        Last major update by: William Wei on 7/12/2024

        """
        # Interpolate EFIT signals to TS time basis
        bminor = interp1(efit_time, bminor, TS_time)
        z0 = interp1(efit_time, z0, TS_time)

        # Calculate Te peaking factor
        Te_PF = np.full(len(TS_time), np.nan)
        (itimes,) = np.where((TS_time > 0) & (TS_time < times[-1]))
        for itime in itimes:
            TS_Te_arr = TS_Te[:, itime]
            (indx,) = np.where(TS_Te_arr > 0)
            if len(indx) < 10:
                continue
            TS_Te_arr = TS_Te_arr[indx]
            TS_z_arr = TS_z[indx]
            sorted_indx = np.argsort(TS_z_arr)
            TS_z_arr = TS_z_arr[sorted_indx]
            TS_Te_arr = TS_Te_arr[sorted_indx]
            # Create equal-spacing array of TS_z_arr and interpolate TS_Te_arr on it
            # Skip if there's no EFIT zmagx data
            if np.isnan(z0[itime]):
                continue
            z_arr_equal_spacing = np.linspace(z0[itime], TS_z_arr[-1], len(TS_z_arr))
            Te_arr_equal_spacing = interp1(TS_z_arr, TS_Te_arr, z_arr_equal_spacing)
            # Calculate Te_PF
            (core_index,) = np.where(
                np.array((z_arr_equal_spacing - z0[itime]) < 0.2 * abs(bminor[itime]))
            )
            if len(core_index) < 2:
                continue
            Te_PF[itime] = np.mean(Te_arr_equal_spacing[core_index]) / np.mean(
                Te_arr_equal_spacing
            )

        # TODO: Calculate ne and pressure peaking factors

        # Interpolate peaking factors to the requested time basis
        Te_PF = interp1(TS_time, Te_PF, times, "linear")
        return {
            "ne_peaking": [np.nan],
            "te_peaking": Te_PF,
            "pressure_peaking": [np.nan],
        }

    @staticmethod
    @physics_method(
        columns=["ne_peaking", "te_peaking", "pressure_peaking"],
        tokamak=Tokamak.CMOD,
    )
    def _get_peaking_factors(params: PhysicsMethodParams):
        nan_output = {
            "ne_peaking": [np.nan],
            "te_peaking": [np.nan],
            "pressure_peaking": [np.nan],
        }
        # Ignore shots on the blacklist
        if CmodPhysicsMethods.is_on_blacklist(params.shot_id):
            return nan_output
        # Fetch data
        try:
            # Get EFIT geometry data
            z0 = 0.01 * params.mds_conn.get_data(
                r"\efit_aeqdsk:zmagx", tree_name="_efit_tree"
            )
            kappa = params.mds_conn.get_data(
                r"\efit_aeqdsk:kappa", tree_name="_efit_tree"
            )
            aminor, efit_time = params.mds_conn.get_data_with_dims(
                r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
            )
            bminor = aminor * kappa
            # Get Te data and TS time basis
            node_ext = ".yag_new.results.profiles"
            TS_Te_core, TS_time = params.mds_conn.get_data_with_dims(
                f"{node_ext}:te_rz", tree_name="electrons"
            )
            # Convert keV to eV
            TS_Te_core = TS_Te_core * 1000
            TS_Te_edge = params.mds_conn.get_data(r"\ts_te")
            # Convert eV to Kelvin
            TS_Te = np.concatenate((TS_Te_core, TS_Te_edge)) * 11600
            TS_z_core = params.mds_conn.get_data(
                f"{node_ext}:z_sorted", tree_name="electrons"
            )
            TS_z_edge = params.mds_conn.get_data(r"\fiber_z", tree_name="electrons")
            TS_z = np.concatenate((TS_z_core, TS_z_edge))
            # Make sure that there are equal numbers of edge position and edge temperature points
            if len(TS_z_edge) != TS_Te_edge.shape[0]:
                params.logger.warning(
                    f"[Shot {params.shot_id}]: TS edge data and z positions are not the same length for shot"
                )
                return nan_output
            # TODO: Get ne data
            TS_ne = np.full(len(params.times), np.nan)
        except mdsExceptions.MdsException as e:
            return nan_output
        return CmodPhysicsMethods.get_peaking_factors(
            params.times, TS_time, TS_Te, TS_ne, TS_z, efit_time, bminor, z0
        )

    @staticmethod
    @physics_method(
        columns=["prad_peaking"],
        tokamak=Tokamak.CMOD,
    )
    def _get_prad_peaking(params: PhysicsMethodParams):
        prad_peaking = np.full(len(params.times), np.nan)
        nan_output = {"prad_peaking": prad_peaking}
        try:
            r0 = 0.01 * params.mds_conn.get_data(
                r"\efit_aeqdsk:rmagx", tree_name="_efit_tree"
            )
            z0 = 0.01 * params.mds_conn.get_data(
                r"\efit_aeqdsk:zmagx", tree_name="_efit_tree"
            )
            aminor, efit_time = params.mds_conn.get_data_with_dims(
                r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
            )
        except mdsExceptions.MdsException as e:
            params.logger.debug(f"[Shot {params.shot_id}]: Failed to get efit data")
            return nan_output
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
        except mdsExceptions.MdsException as e:
            params.logger.debug(f"[Shot {params.shot_id}]: Failed to get AXA data")
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
        except mdsExceptions.MdsException as e:
            params.logger.debug(f"[Shot {params.shot_id}]: Failed to get AXJ data")
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

    @staticmethod
    def get_sxr_parameters():
        pass

    # TODO: get more accurate description of soft x-ray data
    @staticmethod
    @physics_method(columns=["sxr"], tokamak=Tokamak.CMOD)
    def _get_sxr_data(params: PhysicsMethodParams):
        try:
            sxr, t_sxr = params.mds_conn.get_data_with_dims(
                r"\top.brightnesses.array_1:chord_16",
                tree_name="xtomo",
                astype="float64",
            )
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_id}]: Failed to get SXR data returning NaNs"
            )
            params.logger.debug(f"[Shot {params.shot_id}]: {traceback.format_exc()}")
            return {"sxr": [np.nan]}
        sxr = interp1(t_sxr, sxr, params.times)
        return {"sxr": sxr}

    @staticmethod
    def is_on_blacklist(shot_id: int) -> bool:
        """TODO why will these shots cause `_get_peaking_factors`,
        `_get_peaking_factors_no_tci`, and `_get_edge_parameters` to fail?
        """
        if (
            (shot_id > 1120000000 and shot_id < 1120213000)
            or (shot_id > 1140000000 and shot_id < 1140227000)
            or (shot_id > 1150000000 and shot_id < 1150610000)
            or (shot_id > 1160000000 and shot_id < 1160303000)
        ):
            return True
        return False

class CmodTearingMethods:

    @staticmethod
    @physics_method(
        columns=["n_equal_1_std"],
        tokamak=Tokamak.CMOD,
    )
    def _get_n_equal_1_std(params: PhysicsMethodParams):
        """ Calculate n=1 STD from magnetics data.

        This uses the set of Mirnov coils that are not integrated
        (Those that end in k)
        These are distributed in a grid on the outboard side of a limiter.

        At the moment just testing with GH array,
        but there is one on AB as well.
        """

        # Get a single bit of Mirnov coil data
        # Mirnov coils are on a very fast timebase
        # Using the GH signals since that's what Bob did
        # Also have option to use 
        mirnov_names = ["BP02_GHK"]

        # path = r".rf_lim_coils."
        # mirnov_node_names = params.mds_conn.get_data(
        #     path + "nodename", tree_name = "magnetics"
        # )
        # phi = params.mds_conn.get_data(
        #     path + "phi_gh", tree_name = "magnetics"
        # )
        # mirnov_pickup_coeffs = params.mds_conn.get_data(
        #     path + "pickup_coeff", tree_name = "magnetics"
        # )
        # # Only use the coils in mirnov_names
        # _, mirnov_indices, _ = np.intersect1d(
        #     mirnov_node_names, mirnov_names, return_indices=True
        # )
        # mirnov_phi = phi + 360
        # mirnov_pickup_coeffs = mirnov_pickup_coeffs[mirnov_indices]

        # signal = params.mds_conn.get_data(
        #     path + mirnov_names[0], tree_name = "magnetics"
        # )

        path = r"\magnetics::top.active_mhd.signals"
        mirnov_signal, mirnov_times = params.mds_conn.get_data_with_dims(
            path=f"{path}:{mirnov_names[0]}", tree_name="magnetics"
        )

        # Implementation of finding the standard deviation is based on
        # the matlap script get_Mirnov_std.m
        # For each time in the target timebase,
        # Sample the previous 0.001 seconds of the fast mirnov signal
        # and obtain the standard deviation
        # Mirnov is ~5 MHz, so 0.001 seconds has plenty of samples

        target_times = params.times
        time_window = 0.001
        n_equal_1_std = np.zeros(len(target_times))

        for i, time in enumerate(target_times):
            window_indices = np.where(
                (mirnov_times > time - time_window) & (mirnov_times < time)
            )
            # Compute standard deviation ignoring nans
            n_equal_1_std[i] = np.nanstd(mirnov_signal[window_indices])

        return pd.DataFrame(
            {
                "n_equal_1_std": n_equal_1_std,
            }
        )

    @staticmethod
    @physics_method(
        tags=["mirnov_spectrogram"],
        tokamak=Tokamak.CMOD,
    )
    def _get_mirnov_freqs(params: PhysicsMethodParams):
        """Obtain the FFT of the Mirnov coil signals."""

        mirnov_names = ["BP02_GHK"]
        path = r"\magnetics::top.active_mhd.signals"
        mirnov_signal, mirnov_times = params.mds_conn.get_data_with_dims(
            path=f"{path}:{mirnov_names[0]}", tree_name="magnetics"
        )
        
        # Assuming mirnov sample frequency of 5Mhz, and the target timebase is 1 kHz, 
        # we should be able to get 5e6 / 1e3 = 5000 samples per target time
        # This gives us really large frequency bins though,
        # and we're mostly interested in the low frequency range
        # So we'll use a larger number of samples per target time
        mirnov_sample_freq = 5e6 # 5 MHz
        nperseg = 2**15 # 32768 samples per segment

        f, t, Sxx = sp.signal.spectrogram(mirnov_signal, fs=mirnov_sample_freq, nperseg=nperseg)
        # spectrogram thinks time starts at 0
        # so we need to shift the spectrogram timebase to match the mirnov timebase
        t += (mirnov_times[0])

        # Cut the spectrogram to be max 50 kHz
        f_max = 50e3
        f_indices = np.where(f < f_max)
        f = f[f_indices]
        Sxx = Sxx[f_indices]

        # Interpolate the spectrogram onto the target timebase
        # Sxx_interp still has time for columns and frequency for rows
        target_times = params.times
        Sxx_interp = np.zeros((len(f), len(target_times)))
        for i in range(len(f)):
            Sxx_interp[i] = interp1(t, Sxx[i], target_times)

        # Use multi-indexing to store the spectrogram data
        spectrogram_data = pd.DataFrame(Sxx_interp.T, columns=pd.MultiIndex.from_product([["mirnov_spectrogram"], f], names=["", "frequencies"]))

        return spectrogram_data
    
    @staticmethod
    @physics_method(
        tags="sxr_array_1", tokamak=Tokamak.CMOD
    )
    def _get_sxr_chord_data(params: PhysicsMethodParams):
        """ """
        sxr_dataframe = pd.DataFrame()
        for chord in range(0,40):
            try:
                path = f"\\top.brightnesses.array_1:chord_{chord}"
                sxr, t_sxr = params.mds_conn.get_data_with_dims(
                    r"{}".format(path),
                    tree_name="xtomo",
                    astype="float64",
                )
                sxr = interp1(t_sxr, sxr, params.shot_props.times)
            except:
                sxr = np.full(len(params.shot_props.times), np.nan)

            sxr_dataframe[f"sxr_array_1", chord] = sxr
        return sxr_dataframe

    @staticmethod
    @physics_method(
        columns=["sxr_emiss", "sxr_array2_bright", "sxr_array4_bright"], tokamak=Tokamak.CMOD
    )
    def _get_sxr_profile_data(params: PhysicsMethodParams):
        """ """
        try:
            sxr_emiss, t_sxr = params.mds_conn.get_data_with_dims(
                r"\top.results.core:emiss",
                tree_name="xtomo",
                astype="float64",
            )
            sxr_emiss = interp1(t_sxr, sxr_emiss, params.shot_props.times)
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_props.shot_id}]: Failed to get SXR EMISS data returning NaNs"
            )
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
            )
            sxr_emiss = np.full(len(params.shot_props.times), np.nan)

        try:
            sxr_array2_bright, t_sxr = params.mds_conn.get_data_with_dims(
                r"\top.results.array_2:bright",
                tree_name="xtomo",
                astype="float64",
            )
            sxr_array2_bright = interp1(t_sxr, sxr_array2_bright, params.shot_props.times)
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_props.shot_id}]: Failed to get SXR ARRAY 2 BRIGHT data returning NaNs"
            )
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
            )
            sxr_array2_bright = np.full(len(params.shot_props.times), np.nan)

        try:
            sxr_array4_bright, t_sxr = params.mds_conn.get_data_with_dims(
                r"\top.results.array_4:bright",
                tree_name="xtomo",
                astype="float64",
            )
            sxr_array4_bright = interp1(t_sxr, sxr_array4_bright, params.shot_props.times)
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_props.shot_id}]: Failed to get SXR ARRAY 4 BRIGHT data returning NaNs"
            )
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
            )
            sxr_array4_bright = np.full(len(params.shot_props.times), np.nan)
        return pd.DataFrame({"sxr_emiss": sxr_emiss, "sxr_array2_bright": sxr_array2_bright, "sxr_array4_bright": sxr_array4_bright})  


    @staticmethod
    @physics_method(
        tags=["ece_te"], tokamak=Tokamak.CMOD
    )
    def _get_ece_te_data(params: PhysicsMethodParams):
        ece_dataframe = pd.DataFrame()
        for channel_number in range(1, 10):
            try:
                path = f"\\top.ece.gpc_results.te:te{channel_number}"
                ece, t_ece = params.mds_conn.get_data_with_dims(
                    r"{}".format(path),
                    tree_name="electrons",
                    astype="float64",
                )
                ece = interp1(t_ece, ece, params.times)
            except:
                ece = np.full(len(params.times), np.nan)

            ece_dataframe["ece_te", channel_number] = ece

        return ece_dataframe
    
    @staticmethod
    @physics_method(
        tags=["cross_spectrum", "cross_power", "cross_phase", "cross_coherence"], tokamak=Tokamak.CMOD
    )
    def _get_mirnov_spectra(params: PhysicsMethodParams):

        mirnov_names = ["BP02_GHK", "BP02_ABK", "BP14_ABK"]
        path = r"\magnetics::top.active_mhd.signals"
        mirnov_signal_0, mirnov_times = params.mds_conn.get_data_with_dims(
            path=f"{path}:{mirnov_names[0]}", tree_name="magnetics"
        )
        try:
            mirnov_signal_1, mirnov_times = params.mds_conn.get_data_with_dims(
                path=f"{path}:{mirnov_names[1]}", tree_name="magnetics"
            )
        except:
            mirnov_signal_1, mirnov_times = params.mds_conn.get_data_with_dims(
                path=f"{path}:{mirnov_names[2]}", tree_name="magnetics"
            )

        # Assuming mirnov sample frequency of 5Mhz
        mirnov_sample_freq = 5e6

        # Get the complex cross spectrum
        # First, need the short time fourier transform of the signals
        from scipy.signal import ShortTimeFFT
        from scipy.signal.windows import gaussian
        window = gaussian(2**12, std=2**10) # TODO(ZanderKeith), need to make sure sampling window is only in the past
        SFT = ShortTimeFFT(window, hop=2048, fs=mirnov_sample_freq, scale_to="magnitude")
        mirnov_signal_0_fft = SFT.stft(mirnov_signal_0)
        mirnov_signal_1_fft = SFT.stft(mirnov_signal_1)

        freqs = SFT.f
        fft_times = (SFT.delta_t * np.arange(mirnov_signal_0_fft.shape[1])) + mirnov_times[0]
    
        # Cut the fft to be max 50 kHz
        f_max = 50e3
        f_indices = np.where(freqs < f_max)
        freqs = freqs[f_indices]
        mirnov_signal_0_fft = mirnov_signal_0_fft[f_indices]
        mirnov_signal_1_fft = mirnov_signal_1_fft[f_indices]

        # Interpolate the fft's onto the target timebase
        mirnov_signal_0_fft_interp = np.zeros((len(freqs), len(params.times)), dtype=np.complex128)
        mirnov_signal_1_fft_interp = np.zeros((len(freqs), len(params.times)), dtype=np.complex128)
        for i in range(len(freqs)):
            mirnov_signal_0_fft_interp[i] = interp1(fft_times, mirnov_signal_0_fft[i], params.times)
            mirnov_signal_1_fft_interp[i] = interp1(fft_times, mirnov_signal_1_fft[i], params.times)

        # Get the complex cross spectrum
        Cxy = np.zeros((len(freqs), len(params.times)), dtype=np.complex128)
        for i in range(len(freqs)):
                Cxy[i] = mirnov_signal_0_fft_interp[i]*np.conj(mirnov_signal_1_fft_interp[i])

        # Get the cross power
        # TODO(ZanderKeith): There is a dedicated scipy function for cross power, should compare results
        #f, Pxy = sp.signal.csd(mirnov_signal_0, mirnov_signal_1, fs=mirnov_sample_freq, nperseg=nperseg)
        cross_power = np.abs(Cxy)

        # Get the cross phase
        cross_phase = np.angle(Cxy)

        # Get the cross coherence
        cross_coherence = np.zeros((len(freqs), len(params.times)))
        for i in range(len(freqs)):
            a = cross_power[i]**2
            b = (np.abs(mirnov_signal_0_fft_interp[i])**2 * np.abs(mirnov_signal_1_fft_interp[i])**2)
            mask = b != 0
            result = np.true_divide(a, b, where=mask)
            result[~mask] = np.nan
            cross_coherence[i] = result

        # Use multi-indexing to store the cross spectrum data
        cross_spectrum_data = pd.DataFrame(Cxy.T, columns=pd.MultiIndex.from_product([["cross_spectrum"], freqs], names=["", "frequencies"]))
        cross_power_data = pd.DataFrame(cross_power.T, columns=pd.MultiIndex.from_product([["cross_power"], freqs], names=["", "frequencies"]))
        cross_phase_data = pd.DataFrame(cross_phase.T, columns=pd.MultiIndex.from_product([["cross_phase"], freqs], names=["", "frequencies"]))
        cross_coherence_data = pd.DataFrame(cross_coherence.T, columns=pd.MultiIndex.from_product([["cross_coherence"], freqs], names=["", "frequencies"]))

        return pd.concat([cross_spectrum_data, cross_power_data, cross_phase_data, cross_coherence_data], axis=1)

    @staticmethod
    @physics_method(
        tags=["partial_flux"], tokamak=Tokamak.CMOD
    )
    def _get_partial_flux(params: PhysicsMethodParams):
        # TODO: look into calibration stuff. The commented out things don't seem to do anything
        #calibration_path = r"\magnetics::top.flux_partial.btor_pickup"
        #calibration_signal = params.mds_conn.get_data(path=calibration_path, tree_name="magnetics")
        
        partial_flux_names = ["F08", "F13", "F14", "F15", "F20", "F27", "F28", "F29", "LM_CD", "LM_EF", "LM_HJ", "LM_KA"]
        path = r"\magnetics::top.flux_partial.signals"

        partial_flux_dataframe = pd.DataFrame()
        for i, partial_flux_coil in enumerate(partial_flux_names):
            try:
                partial_flux_signal, partial_flux_times = params.mds_conn.get_data_with_dims(
                    path=f"{path}:{partial_flux_coil}", tree_name="magnetics"
                )
                #partial_flux_signal = partial_flux_signal - calibration_signal[i]
                partial_flux_signal_interp = interp1(partial_flux_times, partial_flux_signal, params.times)
            except:
                partial_flux_signal_interp = np.full(len(params.times), np.nan)

            partial_flux_name = "partial_flux_" + partial_flux_coil

            # Interpolate the partial flux onto the target timebase
            partial_flux_dataframe[partial_flux_name] = partial_flux_signal_interp

        return partial_flux_dataframe
