#!/usr/bin/env python3

import logging
import sys
import traceback
import warnings
from importlib import resources

import numpy as np
import pandas as pd
import scipy as sp
from MDSplus import mdsExceptions

import disruption_py.data
from disruption_py.settings.shot_data_request import (
    ShotDataRequestParams,
)
from disruption_py.shots.helpers.method_caching import (
    register_method,
)
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import gaussian_fit, interp1, smooth
from disruption_py.utils.utils import safe_cast

warnings.filterwarnings("error", category=RuntimeWarning)


class CModEfitRequests:

    efit_cols = {
        "beta_n": r"\efit_aeqdsk:betan",
        "beta_p": r"\efit_aeqdsk:betap",
        "kappa": r"\efit_aeqdsk:eout",
        "li": r"\efit_aeqdsk:li",
        "upper_gap": r"\efit_aeqdsk:otop/100",
        "lower_gap": r"\efit_aeqdsk:obott/100",
        "q0": r"\efit_aeqdsk:q0",
        "qstar": r"\efit_aeqdsk:qstar",
        "q95": r"\efit_aeqdsk:q95",
        "v_loop_efit": r"\efit_aeqdsk:vloopt",
        "Wmhd": r"\efit_aeqdsk:wplasm",
        "ssep": r"\efit_aeqdsk:ssep/100",
        "n_over_ncrit": r"-\efit_aeqdsk:xnnc",
        "tritop": r"\efit_aeqdsk:doutu",
        "tribot": r"\efit_aeqdsk:doutl",
        "a_minor": r"\efit_aeqdsk:aminor",
        "rmagx": r"\efit_aeqdsk:rmagx",  # TODO: change units to [m] (current [cm])
        "chisq": r"\efit_aeqdsk:chisq",
    }

    # EFIT column names for data before 2000 TODO: confirm with Bob that these are the right back-ups and make sure that these are similar to standard EFIT columns
    efit_cols_pre_2000 = {
        "a_minor": r"\efit_aeqdsk:aout",
        "li": r"\efit_aeqdsk:ali",
        "q0": r"\efit_aeqdsk:qqmagx",
        "qstar": r"\efit_aeqdsk:qsta",
        "q95": r"\efit_aeqdsk:qsib",  # Not sure about this one
    }

    efit_derivs = {"beta_p": "dbetap_dt", "li": "dli_dt", "Wmhd": "dWmhd_dt"}

    @staticmethod
    @register_method(
        columns=[
            *efit_cols.keys(),
            *efit_cols_pre_2000.keys(),
            *efit_derivs.keys(),
            "V_surf",
            "v_loop_efit",
            "beta_n",
        ],
        tokamak=Tokamak.CMOD,
    )
    def _get_EFIT_parameters(params: ShotDataRequestParams):

        efit_time = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
        )  # [s]
        efit_data = dict()

        # Get data from each of the columns in efit_cols one at a time
        for param in CModEfitRequests.efit_cols:
            try:
                # If shot before 2000 and the param is in efit_cols_pre_2000
                if (
                    params.shot_props.shot_id <= 1000000000
                    and param not in CModEfitRequests.efit_cols_pre_2000.keys()
                ):
                    efit_data[param] = params.mds_conn.get_data(
                        CModEfitRequests.efit_cols_pre_2000[param],
                        tree_name="_efit_tree",
                        astype="float64",
                    )
                else:
                    efit_data[param] = params.mds_conn.get_data(
                        CModEfitRequests.efit_cols[param],
                        tree_name="_efit_tree",
                        astype="float64",
                    )
            except:
                params.logger.warning(
                    f"[Shot {params.shot_props.shot_id}]: Unable to get {param} from EFIT tree"
                )
                params.logger.debug(
                    f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
                )
                efit_data[param] = np.full(len(efit_time), np.nan)
                pass

        for param in CModEfitRequests.efit_derivs:
            efit_data[CModEfitRequests.efit_derivs[param]] = np.gradient(
                efit_data[param], efit_time, edge_order=1
            )

        # Get data for V_surf := deriv(\ANALYSIS::EFIT_SSIBRY)*2*pi
        try:
            ssibry = params.mds_conn.get_data(
                r"\efit_geqdsk:ssibry", tree_name="_efit_tree", astype="float64"
            )
            efit_data["V_surf"] = np.gradient(ssibry, efit_time) * 2 * np.pi
        except:
            print("unable to get V_surf")
            efit_data["V_surf"] = np.full(len(efit_time), np.nan)
            pass

        # For shots before 2000, adjust units of aminor, compute beta_n and v_loop
        if params.shot_props.shot_id <= 1000000000:

            # Adjust aminor units
            efit_data["aminor"] = efit_data["aminor"] / 100  # [cm] to [m]

            # Get data for v_loop --> deriv(\ANALYSIS::EFIT_SSIMAG)*$2pi (not totally sure on this one)
            try:  # TODO: confirm this
                ssimag = params.mds_conn.get_data(
                    r"\efit_geqdsk:ssimag", tree_name="_efit_tree", astype="float64"
                )
                efit_data["v_loop_efit"] = np.gradient(ssimag, efit_time) * 2 * np.pi
            except:
                print("unable to get v_loop_efit")
                efit_data["v_loop_efit"] = np.full(len(efit_time), np.nan)
                pass

            # Compute beta_n
            beta_t = params.mds_conn.get_data(
                r"\efit_aeqdsk:betat", tree_name="_efit_tree", astype="float64"
            )
            efit_data["beta_n"] = np.reciprocal(
                np.reciprocal(beta_t) + np.reciprocal(efit_data["beta_p"])
            )

        if not np.array_equal(params.shot_props.times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(
                    efit_time, efit_data[param], params.shot_props.times
                )

        return pd.DataFrame(efit_data)

    @staticmethod
    def efit_check(params: ShotDataRequestParams):
        """
        # TODO: Get description from Jinxiang
        """
        values = []
        for expr in [
            r"_lf=\analysis::efit_aeqdsk:lflag",
            r"_l0=((sum(_lf,1) - _lf[*,20] - _lf[*,1])==0)",
            r"_n=\analysis::efit_fitout:nitera,(_l0 and (_n>4))",
        ]:
            values.append(params.mds_conn.get(expr, tree_name="analysis"))
        _n = values[2].data()
        valid_indices = np.nonzero(_n)
        (times,) = params.mds_conn.get_dims(
            r"\analysis::efit_aeqdsk:lflag", tree_name="analysis"
        )
        return valid_indices, times[valid_indices]


class BasicCmodRequests:
    @staticmethod
    @register_method(
        populate=False,
        cache_between_threads=False,
        tokamak=Tokamak.CMOD,
    )
    def get_active_wire_segments(params: ShotDataRequestParams):
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
    @register_method(columns=["time_until_disrupt"], tokamak=Tokamak.CMOD)
    def _get_time_until_disrupt(params: ShotDataRequestParams):
        time_until_disrupt = np.full(len(params.shot_props.times), np.nan)
        if params.shot_props.disrupted:
            time_until_disrupt = (
                params.shot_props.disruption_time - params.shot_props.times
            )
        return pd.DataFrame({"time_until_disrupt": time_until_disrupt})

    @staticmethod
    def get_ip_parameters(times, ip, magtime, ip_prog, pcstime):
        """Calculates actual and programmed current as well as their derivatives and difference.

        The time derivatives are useful for discriminating between rampup, flattop, and rampdown.

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
        return pd.DataFrame(
            {
                "ip": ip,
                "dip_dt": dip,
                "dip_smoothed": dip_smoothed,
                "ip_prog": ip_prog,
                "dipprog_dt": dipprog_dt,
                "ip_error": ip_error,
            }
        )

    @staticmethod
    @register_method(
        columns=["ip", "dip_dt", "dip_smoothed", "ip_prog", "dipprog_dt", "ip_error"],
        tokamak=Tokamak.CMOD,
    )
    def _get_ip_parameters(params: ShotDataRequestParams):
        # Automatically generated
        active_segments = BasicCmodRequests.get_active_wire_segments(params=params)

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
                                f"[Shot {params.shot_props.shot_id}]: Error getting PID gains for wire {wire_index}"
                            ]
                        )
                        params.logger.debug(
                            [
                                f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
                            ]
                        )
                    break  # Break out of wire_index loop
        ip, magtime = params.mds_conn.get_data_with_dims(
            r"\ip", tree_name="magnetics", astype="float64"
        )
        return BasicCmodRequests.get_ip_parameters(
            params.shot_props.times, ip, magtime, ip_prog, pcstime
        )

    @staticmethod
    def get_z_parameters(times, z_prog, pcstime, z_error_without_ip, ip, dpcstime):
        """Get values of Z_error, Z_prog, and derived signals from plasma control system (PCS).

        Z_prog is the programmed vertical position of the plasma current centroid, and Z_error is the difference
        between the actual position and that requested (Z_error = Z_cur -
        Z_prog). Thus, the actual (estimated) position, Z_cur, can be calculated.
        And the vertical velocity, v_z, can be taken from the time derivative,
        and the product z_times_v_z ( = Z_cur * v_z) is also calculated.

        Parameters
        ----------
        times : array_like
            Time array for the shot.
        z_prog : array_like
            Programmed vertical position of the plasma current centroid.
        pcstime : array_like
            Time array for the programmed vertical position of the plasma current centroid.
        z_error_without_ip : array_like
            Difference between the actual and programmed vertical position of the plasma current centroid.
        ip : array_like
            Actual plasma current.
        dpcstime : array_like
            Time array for the actual plasma current.

        Returns
        -------
        z_error : array_like
            Difference between the actual and programmed vertical position of the plasma current centroid.
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
        return pd.DataFrame(
            {
                "z_error": z_error,
                "z_prog": z_prog,
                "zcur": z_cur,
                "v_z": v_z,
                "z_times_v_z": z_times_v_z,
            }
        )

    @staticmethod
    @register_method(
        columns=["z_error", "z_prog", "zcur", "v_z", "z_times_v_z"],
        tokamak=Tokamak.CMOD,
    )
    def _get_z_parameters(params: ShotDataRequestParams):
        pcstime = np.array(np.arange(-4, 12.383, 0.001))
        z_prog = np.empty(pcstime.shape)
        z_prog.fill(np.nan)
        z_prog_temp = z_prog.copy()
        z_wire_index = -1
        active_wire_segments = BasicCmodRequests.get_active_wire_segments(params=params)

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
                            f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
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
        # The value of Z_error we read is not in the units we want. It must be *divided* by a factor AND *divided* by the plasma current.
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
        if params.shot_props.shot_id > 1150101000:
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
        return BasicCmodRequests.get_z_parameters(
            params.shot_props.times, z_prog, pcstime, z_error_without_ip, ip, dpcstime
        )

    @staticmethod
    def get_ohmic_parameters(
        times, v_loop, v_loop_time, li, efittime, dip_smoothed, ip
    ):
        """Calculate the ohmic power from the loop voltage, inductive voltage, and plasma current.

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
        R0 = 0.68  # For simplicity, we use R0 = 0.68 m, but we could use \efit_aeqdsk:rmagx
        inductance = 4.0 * np.pi * 1.0e-7 * R0 * li / 2.0
        v_loop = interp1(v_loop_time, v_loop, times)
        inductance = interp1(efittime, inductance, times)
        v_inductive = inductance * dip_smoothed
        v_resistive = v_loop - v_inductive
        p_ohm = ip * v_resistive
        return pd.DataFrame({"p_oh": p_ohm, "v_loop": v_loop})

    @staticmethod
    @register_method(
        columns=["p_oh", "v_loop"],
        tokamak=Tokamak.CMOD,
    )
    def _get_ohmic_parameters(params: ShotDataRequestParams):
        v_loop, v_loop_time = params.mds_conn.get_data_with_dims(
            r"\top.mflux:v0", tree_name="analysis", astype="float64"
        )
        if len(v_loop_time) <= 1:
            return pd.DataFrame(
                {
                    "p_oh": np.zeros(len(params.shot_props.times)),
                    "v_loop": np.zeros(len(params.shot_props.times)),
                }
            )
        li, efittime = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:li", tree_name="_efit_tree", astype="float64"
        )
        ip_parameters = BasicCmodRequests._get_ip_parameters(params=params)
        return BasicCmodRequests.get_ohmic_parameters(
            params.shot_props.times,
            v_loop,
            v_loop_time,
            li,
            efittime,
            ip_parameters["dip_smoothed"],
            ip_parameters["ip"],
        )

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
        return pd.DataFrame(
            {
                "p_rad": p_rad,
                "dprad_dt": dprad,
                "p_lh": p_lh,
                "p_icrf": p_icrf,
                "p_input": p_input,
                "radiated_fraction": rad_fraction,
            }
        )

    @staticmethod
    @register_method(
        columns=["p_rad", "dprad_dt", "p_lh", "p_icrf", "p_input", "radiated_fraction"],
        tokamak=Tokamak.CMOD,
    )
    def _get_power(params: ShotDataRequestParams):
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
        p_oh = BasicCmodRequests._get_ohmic_parameters(params=params)["p_oh"]
        return BasicCmodRequests.get_power(params.shot_props.times, *values, p_oh)

    @staticmethod
    def get_kappa_area(times, aminor, area, a_times):
        return pd.DataFrame(
            {"kappa_area": interp1(a_times, area / (np.pi * aminor**2), times)}
        )

    @staticmethod
    @register_method(columns=["kappa_area"], tokamak=Tokamak.CMOD)
    def _get_kappa_area(params: ShotDataRequestParams):
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
        return BasicCmodRequests.get_kappa_area(
            params.shot_props.times, aminor, area, times
        )

    @staticmethod
    def get_rotation_velocity(times, intensity, time, vel, hirextime):
        """
        Uses spectroscopy graphs of ionized(to hydrogen and helium levels) Argon to calculate velocity. Because of the heat profile of the plasma, suitable measurements are only found near the center
        """
        v_0 = np.empty(len(time))
        # Check that the argon intensity pulse has a minimum count and duration threshold
        valid_indices = np.where(intensity > 1000 & intensity < 10000)
        # Matlab code just multiplies by time delta but that doesn't work in the case where we have different time deltas
        # Instead we sum the time deltas for all valid indices to check the total duration
        if np.sum(time[valid_indices + 1] - time[valid_indices]) >= 0.2:
            v_0 = interp1(hirextime, vel, time)
            # TODO: Determine better threshold
            v_0[np.where(abs(v_0) > 200)] = np.nan
            v_0 *= 1000.0
        v_0 = interp1(time, v_0, times)
        return pd.DataFrame({"v_0": v_0})

    # TODO: Calculate v_mid
    @staticmethod
    @register_method(columns=["v_0"], tokamak=Tokamak.CMOD)
    def _get_rotation_velocity(params: ShotDataRequestParams):
        with resources.path(disruption_py.data, "lock_mode_calib_shots.txt") as fio:
            calibrated = pd.read_csv(fio)
        # Check to see if shot was done on a day where there was a locked
        # mode HIREX calibration by cross checking with list of calibrated
        # runs. If not calibrated, return NaN outputs.
        if params.shot_props.shot_id not in calibrated:
            v_0 = np.empty(len(params.shot_props.times))
            v_0.fill(np.nan)
            return pd.DataFrame({"v_0": v_0})
        try:
            intensity, time = params.mds_conn.get_data_with_dims(
                ".hirex_sr.analysis.a:int", tree_name="spectroscopy", astype="float64"
            )
            vel, hirextime = params.mds_conn.get_data_with_dims(
                ".hirex_sr.analysis.a:vel", tree_name="spectroscopy", astype="float64"
            )
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_props.shot_id}]: Failed to open necessary tress for rotational velocity calculations."
            )
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
            )
            v_0 = np.empty(len(params.shot_props.times))
            v_0.fill(np.nan)
            return pd.DataFrame({"v_0": v_0})
        return BasicCmodRequests.get_rotation_velocity(
            params.shot_props.times, intensity, time, vel, hirextime
        )

    # TODO: Split into static and instance method
    @staticmethod
    def get_n_equal_1_amplitude():
        pass

    # TODO: Try catch failure to get BP13 sensors
    @staticmethod
    @register_method(
        columns=["n_equal_1_mode", "n_equal_1_normalized", "n_equal_1_phase", "BT"],
        tokamak=Tokamak.CMOD,
    )
    def _get_n_equal_1_amplitude(params: ShotDataRequestParams):
        """Calculate n=1 amplitude and phase.

        This method uses the four BP13 Bp sensors near the midplane on the outboard vessel
        wall.  The calculation is done by using a least squares fit to an
        expansion in terms of n = 0 & 1 toroidal harmonics.  The BP13 sensors are
        part of the set used for plasma control and equilibrium reconstruction,
        and their signals have been analog integrated (units: tesla), so they
        don't have to be numerically integrated.  These four sensors were working
        well in 2014, 2015, and 2016.  I looked at our locked mode MGI run on
        1150605, and the different applied A-coil phasings do indeed show up on
        the n=1 signal.

        N=1 toroidal assymmetry in the magnetic fields
        """
        n_equal_1_amplitude = np.empty(len(params.shot_props.times))
        n_equal_1_amplitude.fill(np.nan)
        n_equal_1_normalized = n_equal_1_amplitude.copy()
        n_equal_1_phase = n_equal_1_amplitude.copy()
        # These sensors are placed toroidally around the machine. Letters refer to the 2 ports the sensors were placed between.
        bp13_names = ["BP13_BC", "BP13_DE", "BP13_GH", "BP13_JK"]
        bp13_signals = np.empty((len(params.shot_props.times), len(bp13_names)))

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
        # Toroidal power supply takes time to turn on, from ~ -1.8 and should be on by t=-1. So pick the time before that to calculate baseline
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
                        f"[Shot {params.shot_props.shot_id}] Only one data point for {bp13_names[i]} Returning nans."
                    )
                    return n_equal_1_amplitude, n_equal_1_normalized, n_equal_1_phase
                baseline = np.mean(signal[baseline_indices])
                signal = signal - baseline
                signal = signal - bp13_btor_pickup_coeffs[i] * btor
                bp13_signals[:, i] = interp1(t_mag, signal, params.shot_props.times)
            except mdsExceptions.TreeNODATA as e:
                params.logger.warning(
                    f"[Shot {params.shot_props.shot_id}] No data for {bp13_names[i]}"
                )
                params.logger.debug(f"[Shot {params.shot_props.shot_id}] {e}")
                valid_sensors = False
        # TODO: Examine edge case behavior of sign
        polarity = np.sign(np.mean(btor))
        btor_magnitude = btor * polarity
        btor_magnitude = interp1(t_mag, btor_magnitude, params.shot_props.times)
        btor = interp1(t_mag, btor, params.shot_props.times)  # Interpolate BT with sign
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
        return pd.DataFrame(
            {
                "n_equal_1_mode": n_equal_1_amplitude,
                "n_equal_1_normalized": n_equal_1_normalized,
                "n_equal_1_phase": n_equal_1_phase,
                "BT": btor,
            }
        )

    @staticmethod
    def get_densities(times, n_e, t_n, ip, t_ip, a_minor, t_a):
        if len(n_e) != len(t_n):
            nan_arr = np.empty(len(times))
            nan_arr.fill(np.nan)
            return pd.DataFrame(
                {
                    "n_e": nan_arr,
                    "dn_dt": nan_arr.copy(),
                    "Greenwald_fraction": nan_arr.copy(),
                }
            )
        # get the gradient of n_E
        dn_dt = np.gradient(n_e, t_n)
        n_e = interp1(t_n, n_e, times)
        dn_dt = interp1(t_n, dn_dt, times)
        ip = -ip / 1e6  # Convert from A to MA and take positive value
        ip = interp1(t_ip, ip, times)
        a_minor = interp1(t_a, a_minor, times)
        # make sure aminor is not 0 or less than 0
        a_minor[a_minor <= 0] = 0.001
        n_G = ip / (np.pi * a_minor**2) * 1e20  # Greenwald density in m ^-3
        g_f = abs(n_e / n_G)
        return pd.DataFrame({"n_e": n_e, "dn_dt": dn_dt, "Greenwald_fraction": g_f})

    @staticmethod
    @register_method(
        columns=["n_e", "dn_dt", "Greenwald_fraction"],
        tokamak=Tokamak.CMOD,
    )
    def _get_densities(params: ShotDataRequestParams):
        try:
            # Line-integrated density
            n_e, t_n = params.mds_conn.get_data_with_dims(
                r".tci.results:nl_04", tree_name="electrons", astype="float64"
            )
            # Divide by chord length of ~0.6m to get line averaged density.
            # For future refernce, chord length is stored in .01*\analysis::efit_aeqdsk:rco2v[3,*]
            n_e = np.squeeze(n_e) / 0.6
            ip, t_ip = params.mds_conn.get_data_with_dims(
                r"\ip", tree_name="magnetics", astype="float64"
            )
            a_minor, t_a = params.mds_conn.get_data_with_dims(
                r"\efit_aeqdsk:aminor", tree_name="analysis", astype="float64"
            )
        except Exception as e:
            params.logger.debug(f"[Shot {params.shot_props.shot_id}] {e}")
            params.logger.warning(f"[Shot {params.shot_props.shot_id}] No density data")
            # TODO: Handle this case
            raise NotImplementedError(
                "Can't currently handle failure of grabbing density data"
            )
        return BasicCmodRequests.get_densities(
            params.shot_props.times, n_e, t_n, ip, t_ip, a_minor, t_a
        )

    @staticmethod
    def get_efc_current(times, iefc, t_iefc):
        return pd.DataFrame({"I_efc": interp1(t_iefc, iefc, times, "linear")})

    @staticmethod
    @register_method(columns=["I_efc"], tokamak=Tokamak.CMOD)
    def _get_efc_current(params: ShotDataRequestParams):
        try:
            iefc, t_iefc = params.mds_conn.get_data_with_dims(
                r"\efc:u_bus_r_cur", tree_name="engineering"
            )
        except Exception as e:
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}] {traceback.format_exc()}"
            )
            return pd.DataFrame({"I_efc": np.empty(len(params.shot_props.times))})
        return BasicCmodRequests.get_efc_current(params.shot_props.times, iefc, t_iefc)

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
        return pd.DataFrame({"Te_width": te_hwm})

    @staticmethod
    @register_method(columns=["Te_width"], tokamak=Tokamak.CMOD)
    def _get_Ts_parameters(params: ShotDataRequestParams):
        # TODO: Guassian vs parabolic fit for te profile
        te_hwm = np.empty((len(params.shot_props.times)))

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
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}] {traceback.format_exc()}"
            )
            te_hwm.fill(np.nan)
            return pd.DataFrame({"Te_width": te_hwm})
        return BasicCmodRequests.get_Ts_parameters(
            params.shot_props.times, ts_data, ts_time, ts_z
        )

    # TODO: Finish
    @staticmethod
    def get_peaking_factors(times, TS_time, ne_PF, Te_PF, pressure_PF):
        # ne_PF = interp1(TS_time, ne_PF, times, 'linear')
        # Te_PF = interp1(TS_time, Te_PF, times, 'linear')
        # pressure_PF = interp1(TS_time, pressure_PF, times, 'linear')
        pass

    @staticmethod
    @register_method(
        columns=["ne_peaking", "Te_peaking", "pressure_peaking"],
        tokamak=Tokamak.CMOD,
    )
    def _get_peaking_factors(params: ShotDataRequestParams):
        ne_PF = np.full(len(params.shot_props.times), np.nan)
        Te_PF = ne_PF.copy()
        pressure_PF = ne_PF.copy()
        if (
            (
                params.shot_props.shot_id > 1120000000
                and params.shot_props.shot_id < 1120213000
            )
            or (
                params.shot_props.shot_id > 1140000000
                and params.shot_props.shot_id < 1140227000
            )
            or (
                params.shot_props.shot_id > 1150000000
                and params.shot_props.shot_id < 1150610000
            )
            or (
                params.shot_props.shot_id > 1160000000
                and params.shot_props.shot_id < 1160303000
            )
        ):
            # Ignore shots on the blacklist
            return pd.DataFrame(
                {
                    "ne_peaking": ne_PF,
                    "Te_peaking": Te_PF,
                    "pressure_peaking": pressure_PF,
                }
            )
        try:
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
            node_ext = ".yag_new.results.profiles"
            # nl_ts1, nl_ts2, nl_tci1, nl_tci2, _, _ = ThomsonDensityMeasure.compare_ts_tci(params, nlnum=4)
            TS_te, TS_time = params.mds_conn.get_data_with_dims(
                f"{node_ext}:te_rz", tree_name="electrons"
            )
            TS_te = TS_te * 1000 * 11600
            tets_edge = params.mds_conn.get_data(r"\ts_te") * 11600
            TS_te = np.concatenate((TS_te, tets_edge))
            TS_z = params.mds_conn.get_data(
                f"{node_ext}:z_sorted", tree_name="electrons"
            )
            zts_edge = params.mds_conn.get_data(r"\fiber_z", tree_name="electrons")
            TS_z = np.concatenate((TS_z, zts_edge))
            if len(zts_edge) != tets_edge.shape[1]:
                return pd.DataFrame(
                    {
                        "ne_peaking": ne_PF,
                        "Te_peaking": Te_PF,
                        "pressure_peaking": pressure_PF,
                    }
                )
            Te_PF = Te_PF[: len(TS_time)]
            itimes = np.where(TS_time > 0 & TS_time < params.shot_props.times[-1])
            bminor = interp1(efit_time, bminor, TS_time)
            z0 = interp1(efit_time, z0, TS_time)
            for i in range(len(itimes)):
                Te_arr = TS_te[itimes[i], :]
                indx = np.where(Te_arr > 0)
                if len(indx) < 10:
                    continue
                Te_arr = Te_arr[indx]
                TS_z_arr = TS_z[indx]
                sorted_indx = np.argsort(TS_z_arr)
                Ts_z_arr = Ts_z_arr[sorted_indx]
                Te_arr = Te_arr[sorted_indx]
                z_arr = np.linspace(z0[itimes[i]], TS_z_arr[-1], len(Ts_z_arr))
                Te_arr = interp1(TS_z_arr, Te_arr, z_arr)
                core_index = np.where(
                    z_arr
                    < (z0[itimes[i]] + 0.2 * bminor[itimes[i]]) & z_arr
                    > (z0[itimes[i]] - 0.2 * bminor[itimes[i]])
                )
                if len(core_index) < 2:
                    continue
                Te_PF[itimes[i]] = np.mean(Te_arr[core_index]) / np.mean(Te_arr)
            Te_PF = interp1(TS_time, Te_PF, params.shot_props.times)
            calib = np.nan
            # TODO(lajz): fix
            return BasicCmodRequests.get_Ts_parameters(
                params.shot_props.times, TS_time, ne_PF, Te_PF, pressure_PF
            )
        except mdsExceptions.MdsException as e:
            return pd.DataFrame(
                {
                    "ne_peaking": ne_PF,
                    "Te_peaking": Te_PF,
                    "pressure_peaking": pressure_PF,
                }
            )

    @staticmethod
    @register_method(
        columns=["prad_peaking"],
        tokamak=Tokamak.CMOD,
    )
    def _get_prad_peaking(params: ShotDataRequestParams):
        prad_peaking = np.full(len(params.shot_props.times), np.nan)
        try:
            # TODO: why use CMOD here?
            r0 = 0.01 * params.mds_conn.get_data(
                r"\efit_aeqdsk:rmagx", tree_name="cmod"
            )
            z0 = 0.01 * params.mds_conn.get_data(
                r"\efit_aeqdsk:zmagx", tree_name="cmod"
            )
            aminor, efit_time = params.mds_conn.get_data_with_dims(
                r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
            )
        except mdsExceptions.MdsException as e:
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: Failed to get efit data"
            )
            return pd.DataFrame({"prad_peaking": prad_peaking})
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
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: Failed to get AXA data"
            )
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
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: Failed to get AXJ data"
            )
        if not (got_axa or got_axj):
            return pd.DataFrame({"prad_peaking": prad_peaking})
        a_minor = interp1(efit_time, aminor, params.shot_props.times)
        r0 = interp1(efit_time, r0, params.shot_props.times)
        z0 = interp1(efit_time, z0, params.shot_props.times)
        if got_axa:
            good_axa = np.where(good_axa > 0)[0]
            bright_axa = bright_axa[:, good_axa]
            axa_interp = np.full(
                (bright_axa.shape[1], len(params.shot_props.times)), np.nan
            )
            r_axa = r_axa[good_axa]
            for i in range(bright_axa.shape[1]):
                interped = interp1(t_axa.T, bright_axa[:, i], params.shot_props.times.T)
                indx = np.where(interped < 0)
                interped[indx] = np.nan
                axa_interp[i, :] = interped
        if got_axj:
            good_axj = np.where(good_axj > 0)[0]
            bright_axj = bright_axj[:, good_axj]
            axj_interp = np.full(
                (bright_axj.shape[1], len(params.shot_props.times)), np.nan
            )
            r_axj = r_axj[good_axj]
            for i in range(bright_axj.shape[1]):
                interped = interp1(t_axj.T, bright_axj[:, i], params.shot_props.times.T)
                indx = np.where(interped < 0)
                interped[indx] = np.nan
                axj_interp[i, :] = interped
        for i in range(len(params.shot_props.times)):
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
            try:
                prad_peaking[i] = np.nanmean(core_radiation) / np.nanmean(all_radiation)
            except:
                prad_peaking[i] = np.nan
        return pd.DataFrame({"prad_peaking": prad_peaking})

    @staticmethod
    @register_method(
        columns=["ne_peaking", "Te_peaking", "pressure_peaking"],
        tags=["experimental"],
        tokamak=Tokamak.CMOD,
    )
    def _get_peaking_factors_no_tci(params: ShotDataRequestParams):
        # Initialize PFs as empty arrarys
        ne_PF = np.full(len(params.shot_props.times), np.nan)
        Te_PF = ne_PF.copy()
        pressure_PF = ne_PF.copy()
        # Ignore shots on the blacklist
        if (
            (
                params.shot_props.shot_id > 1120000000
                and params.shot_props.shot_id < 1120213000
            )
            or (
                params.shot_props.shot_id > 1140000000
                and params.shot_props.shot_id < 1140227000
            )
            or (
                params.shot_props.shot_id > 1150000000
                and params.shot_props.shot_id < 1150610000
            )
            or (
                params.shot_props.shot_id > 1160000000
                and params.shot_props.shot_id < 1160303000
            )
        ):
            return pd.DataFrame(
                {
                    "ne_peaking": ne_PF,
                    "Te_peaking": Te_PF,
                    "pressure_peaking": pressure_PF,
                }
            )
        try:
            # Get shaping params
            z0 = 0.01 * params.mds_conn.get_data(
                r"\efit_aeqdsk:zmagx", tree_name="cmod"
            )
            kappa = params.mds_conn.get_data(r"\efit_aeqdsk:kappa", tree_name="cmod")
            aminor, efit_time = params.mds_conn.get_data_with_dims(
                r"\efit_aeqdsk:aminor", tree_name="_efit_tree"
            )
            bminor = aminor * kappa  # length of major axis of plasma x-section

            # Get data from TS
            node_ext = ".yag_new.results.profiles"
            # nl_ts1, nl_ts2, nl_tci1, nl_tci2, _, _ = ThomsonDensityMeasure.compare_ts_tci(params, nlnum=4)
            Te_core, Te_time = params.mds_conn.get_data_with_dims(
                f"{node_ext}:te_rz", tree_name="electrons"
            )
            Te_core = Te_core * 1000 * 11600  # Get core TS data
            Te_edge = (
                params.mds_conn.get_data(r"\ts_te", tree_name="electrons") * 11600
            )  # Get edge TS data
            # Concat core and edge data
            Te = np.concatenate((Te_core, Te_edge))
            z_core = params.mds_conn.get_data(
                f"{node_ext}:z_sorted", tree_name="electrons"
            )  # Get z position of core TS points
            # Get z position of edge TS points
            z_edge = params.mds_conn.get_data(r"\fiber_z", tree_name="electrons")
            z = np.concatenate((z_core, z_edge))  # Concat core and edge data
            # Make sure that there are equal numbers of edge position and edge temperature points
            if len(z_edge) != Te_edge.shape[0]:
                params.logger.warning(
                    f"[Shot {params.shot_props.shot_id}]: TS edge data and z positions are not the same length for shot"
                )
                return pd.DataFrame(
                    {
                        "ne_peaking": ne_PF,
                        "Te_peaking": Te_PF,
                        "pressure_peaking": pressure_PF,
                    }
                )
            Te_PF = Te_PF[: len(Te_time)]  # Reshape Te_PF to length of Te_time
            itimes = np.where((Te_time > 0) & (Te_time < params.shot_props.times[-1]))
            node_path = ".yag_new.results.profiles"
            (TS_time,) = params.mds_conn.get_dims(
                node_path + ":te_rz", tree_name="electrons"
            )
            # Interpolate bminor onto desired timebase
            bminor = interp1(efit_time, bminor, TS_time)
            # Interpolate z0 onto desired timebase
            z0 = interp1(efit_time, z0, TS_time)
            for i in range(len(itimes)):
                Te_arr = Te[itimes[i], :]
                indx = np.where(Te_arr > 0)
                if len(indx) < 10:
                    continue
                Te_arr = Te_arr[indx]
                TS_z_arr = z[indx]
                sorted_indx = np.argsort(TS_z_arr)  # Sort by z
                Ts_z_arr = Ts_z_arr[sorted_indx]
                Te_arr = Te_arr[sorted_indx]  # Sort by z
                z_arr = np.linspace(z0[itimes[i]], TS_z_arr[-1], len(Ts_z_arr))
                Te_arr = interp1(TS_z_arr, Te_arr, z_arr)
                core_index = np.where(
                    z_arr
                    < (z0[itimes[i]] + 0.2 * bminor[itimes[i]]) & z_arr
                    > (z0[itimes[i]] - 0.2 * bminor[itimes[i]])
                )
                if len(core_index) < 2:
                    continue
                Te_PF[itimes[i]] = np.mean(Te_arr[core_index]) / np.mean(Te_arr)
            Te_PF = interp1(TS_time, Te_PF, params.shot_props.times)
            calib = np.nan
            # TODO(lajz): fix
            return BasicCmodRequests.get_Ts_parameters(
                params.shot_props.times, TS_time, ne_PF, Te_PF, pressure_PF
            )
        except mdsExceptions.MdsException as e:
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]:{traceback.format_exc()}"
            )
            return pd.DataFrame(
                {
                    "ne_peaking": ne_PF,
                    "Te_peaking": Te_PF,
                    "pressure_peaking": pressure_PF,
                }
            )

    @staticmethod
    def get_sxr_parameters():
        pass

    # TODO: get more accurate description of soft x-ray data
    @staticmethod
    @register_method(columns=["sxr"], tokamak=Tokamak.CMOD)
    def _get_sxr_data(params: ShotDataRequestParams):
        """ """
        sxr = np.full(len(params.shot_props.times), np.nan)
        try:
            sxr, t_sxr = params.mds_conn.get_data_with_dims(
                r"\top.brightnesses.array_1:chord_16",
                tree_name="xtomo",
                astype="float64",
            )
            sxr = interp1(t_sxr, sxr, params.shot_props.times)
        except mdsExceptions.TreeFOPENR as e:
            params.logger.warning(
                f"[Shot {params.shot_props.shot_id}]: Failed to get SXR data returning NaNs"
            )
            params.logger.debug(
                f"[Shot {params.shot_props.shot_id}]: {traceback.format_exc()}"
            )
        return pd.DataFrame({"sxr": sxr})

    @staticmethod
    def get_edge_parameters(times, p_Te, p_ne, edge_rho_min=0.85, edge_rho_max=0.95):
        """Compute the edge Temperature and edge Density signal from the TS.

        Parameters
        ----------
        times : array_like
            The times at which to calculate the edge parameters.
        p_Te : BivariatePlasmaProfile
            The Te measurements [keV] in terms of the time and rho of the measurment.
        p_ne : BivariatePlasmaProfile
            The ne measurements [keV] in terms of the time and rho of the measurment.
        edge_rho_min : float [0,1]
            The rho that defines the minimum of the "edge" region
        edge_rho_max : float [0,1]
            The rho that defines the maximum of the "edge" region

        Returns
        -------
        Te_edge : array_like
            The edge temperature (averaged over the edge region) on the requested timebase.
        ne_edge : array_like
            The edge density (averaged over the edge region) on the requested timebase.

        Original Authors
        ----------------
        Andrew Maris (maris@mit.edu)


        """
        # Base of rho to interpolate onto
        rhobase = np.arange(0, 1, 0.001)

        # Linear interpolate on time and rho
        Te_interpolator = sp.interpolate.LinearNDInterpolator(
            (p_Te.X[:, 0], p_Te.X[:, 1]), p_Te.y
        )
        ne_interpolator = sp.interpolate.LinearNDInterpolator(
            (p_ne.X[:, 0], p_ne.X[:, 1]), p_ne.y
        )

        # Create mesh to compute interpolation over
        timebase_mesh, rhobase_mesh = np.meshgrid(times, rhobase)
        # Compute interpolated values
        Te_interp = Te_interpolator(timebase_mesh, rhobase_mesh)
        ne_interp = ne_interpolator(timebase_mesh, rhobase_mesh)

        plotting = False
        if plotting:
            import matplotlib.pyplot as plt

            plt.ion()
            plt.pcolormesh(timebase_mesh, rhobase_mesh, Te_interp, shading="auto")
            plt.plot(p_Te.X[:, 0], p_Te.X[:, 1], "ok", label="input point")
            plt.legend()
            plt.colorbar()
            plt.show(block=True)

        # Compute Te_edge
        # Make mask for rho in edge region
        rhobase_mesh_mask = (rhobase_mesh >= edge_rho_min) & (
            rhobase_mesh <= edge_rho_max
        )

        # Assert that rho values are indeed in desired range
        assert np.all(rhobase_mesh[rhobase_mesh_mask] >= edge_rho_min)
        assert np.all(rhobase_mesh[rhobase_mesh_mask] <= edge_rho_max)

        # Use mask to get only edge values
        Te_interp_edge = np.where(rhobase_mesh_mask, Te_interp, np.nan)
        ne_interp_edge = np.where(rhobase_mesh_mask, ne_interp, np.nan)

        # Compute edge quantities
        with warnings.catch_warnings():  # Carch warning about taking nanmean of an empty array. This is ok because we want it to return nan for empty arrays
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            Te_edge = np.nanmean(Te_interp_edge, axis=0)
            ne_edge = np.nanmean(ne_interp_edge, axis=0)
        assert len(Te_edge) == len(times)
        assert len(ne_edge) == len(times)

        if plotting:
            plt.plot(times, Te_edge, label="Te_edge")
            plt.plot(times, ne_edge, label="ne_edge")
            plt.legend()
            plt.show(block=True)

        return pd.DataFrame({"Te_edge": Te_edge, "ne_edge": ne_edge})

    @staticmethod
    @register_method(
        tags=["experimental"],
        columns=["Te_edge", "ne_edge"],
        tokamak=Tokamak.CMOD,
    )
    def _get_edge_parameters(params: ShotDataRequestParams):

        try:
            # sys.path.append("/home/sciortino/usr/python3modules/eqtools3")
            sys.path.append("/home/sciortino/usr/python3modules/profiletools3")
            sys.path.append("/home/sciortino/usr/python3modules/gptools3")
            # import eqtools
            import profiletools
        except Exception as e:
            logging.warning("Could not import profiletools or eqtools")
            logging.debug(traceback.format_exc())
            pass

        # Ignore shots on the blacklist
        if (
            (
                params.shot_props.shot_id > 1120000000
                and params.shot_props.shot_id < 1120213000
            )
            or (
                params.shot_props.shot_id > 1140000000
                and params.shot_props.shot_id < 1140227000
            )
            or (
                params.shot_props.shot_id > 1150000000
                and params.shot_props.shot_id < 1150610000
            )
            or (
                params.shot_props.shot_id > 1160000000
                and params.shot_props.shot_id < 1160303000
            )
        ):
            return pd.DataFrame(
                {
                    "Te_edge": np.full(len(params.shot_props.times), np.nan),
                    "ne_edge": np.full(len(params.shot_props.times), np.nan),
                }
            )

        # Range of rho to interpolate over
        rhobase = np.arange(0, 1, 0.001)
        # Get mina and max time from TS tree
        node_path = ".yag_new.results.profiles"
        try:
            (ts_time,) = params.mds_conn.get_dims(
                node_path + ":te_rz", tree_name="electrons"
            )
        except:
            return pd.DataFrame(
                {
                    "Te_edge": np.full(len(params.shot_props.times), np.nan),
                    "ne_edge": np.full(len(params.shot_props.times), np.nan),
                }
            )

        t_min = np.max([0.1, np.min(ts_time)])
        t_max = np.max(ts_time)

        # Get core and edge Thomson profiles over rho := sqrtpsinorm
        p_Te = profiletools.Te(
            params.shot_props.shot_id,
            include=["CTS", "ETS"],
            abscissa="sqrtpsinorm",
            t_min=t_min,
            t_max=t_max,
            remove_zeros=True,
        )
        p_ne = profiletools.ne(
            params.shot_props.shot_id,
            include=["CTS", "ETS"],
            abscissa="sqrtpsinorm",
            t_min=t_min,
            t_max=t_max,
            remove_zeros=True,
        )

        # try:
        #    equal_R = p_ne.X[:,1] == p_Te.X[:,1]
        #    assert np.sum(equal_R) == len(p_ne.X[:,1])
        # except:
        #    raise ValueError('Edge Thomson rhobase differs between ne and Te')
        #    return None, None

        # consider only flux surface on which points were measured, regardless of LFS or HFS
        p_Te.X = np.abs(p_Te.X)
        p_ne.X = np.abs(p_ne.X)

        # set some minimum uncertainties. Recall that units in objects are 1e20m^{-3} and keV
        p_ne.y[p_ne.y <= 0.0] = 0.01  # 10^18 m^-3
        p_Te.y[p_Te.y <= 0.01] = 0.01  # 10 eV
        p_ne.err_y[p_ne.err_y <= 0.01] = 0.01  # 10^18 m^-3
        p_Te.err_y[p_Te.err_y <= 0.02] = 0.02  # 20 eV

        # points in the pedestal that have x uncertainties larger than 0.1 don't help at all
        # do this filtering here because filtering of err_X only works before time-averaging
        p_ne.remove_points(np.logical_and(p_ne.X[:, 1] >= 0.85, p_ne.err_X[:, 1] > 0.1))
        p_Te.remove_points(np.logical_and(p_Te.X[:, 1] >= 0.85, p_Te.err_X[:, 1] > 0.1))

        # cleanup of low Te values
        # TS Te should be >15 eV inside near SOL
        p_Te.remove_points(np.logical_and(p_Te.X[:, 0] < 1.03, p_Te.y < 0.015))

        return BasicCmodRequests.get_edge_parameters(
            params.shot_props.times, p_Te, p_ne
        )

    @staticmethod
    def get_H98():
        pass

    # TODO: Finish
    @staticmethod
    @register_method(
        tags=["experimental"],
        columns=["H98", "Wmhd", "btor", "dWmhd_dt", "p_input"],
        tokamak=Tokamak.CMOD,
    )
    def _get_H98(params: ShotDataRequestParams):
        """Prepare to compute H98 by getting tau_E

        Scaling from eq. 20, ITER Physics Basis Chapter 2 https://iopscience.iop.org/article/10.1088/0029-5515/39/12/302/pdf
        (in s, MA, T, MW, 10^19 m^−3, AMU, m)
        Original Authors
        ----------------
        Andrew Maris (maris@mit.edu)

        """

        # Get parameters for calculating confinement time
        powers_df = BasicCmodRequests._get_power(params=params)
        efit_df = CModEfitRequests._get_EFIT_parameters(params=params)
        density_df = BasicCmodRequests._get_densities(params=params)
        ip_df = BasicCmodRequests._get_ip_parameters(params=params)

        # Get BT

        btor, t_mag = params.mds_conn.get_data_with_dims(
            r"\btor", tree_name="magnetics"
        )  # tmag: [s]
        # Toroidal power supply takes time to turn on, from ~ -1.8 and should be on by t=-1. So pick the time before that to calculate baseline
        baseline_indices = np.where(t_mag <= -1.8)
        btor = btor - np.mean(btor[baseline_indices])
        btor = np.abs(interp1(t_mag, btor, params.shot_props.times))

        ip = np.abs(ip_df.ip) / 1.0e6  # [A] -> [MA]
        n_e = density_df.n_e / 1.0e19  # [m^-3] -> [10^19 m^-3]
        p_input = powers_df.p_input / 1.0e6  # [W] -> [MW]
        dWmhd_dt = efit_df.dWmhd_dt / 1.0e6  # [W] -> [MW]
        Wmhd = efit_df.Wmhd / 1.0e6  # [J] -> [MJ]
        R0 = efit_df.rmagx / 100  # [cm] -> [m]
        # Estimate confinement time
        tau = Wmhd / (p_input - dWmhd_dt)

        # Compute 1998 tau_E scaling, taking A (atomic mass) = 2
        tau_98 = (
            0.0562
            * (n_e**0.41)
            * (2**0.19)
            * (ip**0.93)
            * (R0**1.39)
            * (efit_df.a_minor**0.58)
            * (efit_df.kappa**0.78)
            * (btor**0.15)
            * (p_input**-0.69)
        )
        H98 = tau / tau_98

        return pd.DataFrame(
            {
                "H98": H98,
                "Wmhd": Wmhd,
                "btor": btor,
                "dWmhd_dt": dWmhd_dt,
                "p_input": p_input,
            }
        )


# helper class holding functions for thomson density measures
class ThomsonDensityMeasure:

    # The following methods are translated from IDL code.
    @staticmethod
    def compare_ts_tci(params: ShotDataRequestParams, nlnum=4):
        """
        Comparison between chord integrated Thomson electron density and TCI results.
        """
        core_mult = 1.0
        edge_mult = 1.0
        nl_ts1 = [1e32]
        nl_ts2 = [1e32]
        nl_tci1 = [1e32]
        nl_tci2 = [1e32]
        ts_time1 = [1e32]
        ts_time2 = [1e32]
        (tci_time,) = params.mds_conn.get_dims(
            ".YAG_NEW.RESULTS.PROFILES:NE_RZ", tree_name="electrons"
        )
        tci, tci_t = params.mds_conn.get_data_with_dims(
            f".TCI.RESULTS:NL_{nlnum:02d}", tree_name="electrons"
        )
        nlts, nlts_t = ThomsonDensityMeasure.integrate_ts_tci(nlnum)
        t0 = np.amin(nlts_t)
        t1 = np.amax(nlts_t)
        nyag1, nyag2, indices1, indices2 = ThomsonDensityMeasure.parse_yags(params)
        if nyag1 > 0:
            indices1 += 1
            ts_time1 = tci_time[indices1]
            valid_indices = np.where(ts_time1 >= t0 & ts_time1 <= t1)
            if valid_indices.size > 0:
                nl_tci1 = interp1(tci_t, tci, ts_time1[valid_indices])
                nl_ts1 = interp1(nlts_t, nlts, ts_time1[valid_indices])
                time1 = ts_time1[valid_indices]
        else:
            time1 = -1
        if nyag2 > 0:
            indices2 += 1
            ts_time2 = tci_time[indices2]
            valid_indices = np.where(ts_time2 >= t0 & ts_time2 <= t1)
            if valid_indices.size > 0:
                nl_tci1 = interp1(tci_t, tci, ts_time2[valid_indices])
                nl_ts1 = interp1(nlts_t, nlts, ts_time2[valid_indices])
                time2 = ts_time2[valid_indices]
        else:
            time2 = -1
        return nl_ts1, nl_ts2, nl_tci1, nl_tci2, time1, time2

    @staticmethod
    def parse_yags(params: ShotDataRequestParams):
        nyag1 = params.mds_conn.get_data(r"\knobs:pulses_q", tree_name="electrons")
        nyag2 = params.mds_conn.get_data(r"\knobs:pulses_q_2", tree_name="electrons")
        indices1 = -1
        indices2 = -1
        dark = params.mds_conn.get_data(r"\n_dark_prior", tree_name="electrons")
        ntotal = params.mds_conn.get_data(r"\n_total", tree_name="electrons")
        nt = ntotal - dark
        if nyag1 == 0:
            if nyag2 != 0:
                indices2 = np.arange(nyag2)
        else:
            if nyag2 == 0:
                indices1 = np.arange(nyag1)
            else:
                if nyag1 == nyag2:
                    indices1 = 2 * np.arange(nyag1)
                    indices2 = indices1 + 1
                else:
                    if nyag1 == nyag2:
                        indices1 = 2 * np.arange(nyag1)
                        indices2 = indices1 + 1
                    else:
                        indices1 = 2 * np.arange(nyag1) + (nyag1 > nyag2)
                        indices2 = np.concatenate(
                            (
                                2 * np.arange(nyag2) + (nyag1 < nyag2),
                                2 * nyag2 + np.arange(nyag1 - nyag2 - 1),
                            )
                        )
        v_ind1 = np.where(indices1 < nt)
        if nyag1 > 0 and v_ind1.size > 0:
            indices1 = indices1[v_ind1]
        else:
            indices1 = -1
        v_ind2 = np.where(indices2 < nt)
        if nyag2 > 0 and v_ind2.size > 0:
            indices2 = indices2[v_ind2]
        else:
            indices2 = -1
        return nyag1, nyag2, indices1, indices2

    @staticmethod
    def integrate_ts_tci(params: ShotDataRequestParams, nlnum):
        """
        Integrate Thomson electron density measurement to the line integrated electron density for comparison with two color interferometer (TCI) measurement results
        """
        core_mult = 1.0
        edge_mult = 1.0
        nlts = 1e32
        nlts_t = 1e32
        t, z, n_e, n_e_sig = ThomsonDensityMeasure.map_ts2tci(params, nlnum)
        if z[0, 0] == 1e32:
            return None, None  # TODO: Log and maybe return nan arrs
        nts = len(t)
        nlts_t = t
        nlts = np.full(t.shape, np.nan)
        for i in range(len(nts)):
            ind = np.where(
                np.abs(z[i, :])
                < 0.5 & n_e[i, :]
                > 0 & n_e[i, :]
                < 1e21 & n_e[i, :] / n_e_sig[i, :]
                > 2
            )
            if len(ind) < 3:
                nlts[i] = 0
            else:
                x = z[i, ind]
                y = n_e[i, ind]
                values_uniq, ind_uniq = np.unique(x, return_index=True)
                y = y[ind_uniq]
                nlts[i] = np.trapz(y, x)
        return nlts, nlts_t

    @staticmethod
    def map_ts2tci(params: ShotDataRequestParams, nlnum):
        core_mult = 1.0
        edge_mult = 1.0
        t = [1e32]
        z = [1e32]
        n_e = [1e32]
        n_e_sig = [1e32]
        flag = 1
        valid_indices, efit_times = CModEfitRequests.efit_check()
        ip = params.mds_conn.get_data(r"\ip", "cmod")
        if np.mean(ip) > 0:
            flag = 0
        efit_times = params.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
        )
        t1 = np.amin(efit_times)
        t2 = np.amax(efit_times)
        psia, psia_t = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:SIBDRY", tree_name="_efit_tree"
        )
        psi_0 = params.mds_conn.get(r"\efit_aeqdsk:SIMAGX", tree_name="_efit_tree")
        nets_core, nets_core_t = params.mds_conn.get_data_with_dims(
            ".YAG_NEW.RESULTS.PROFILES:NE_RZ", tree_name="electrons"
        )
        nets_core_err = params.mds_conn.get_data(
            ".YAG_NEW.RESULTS.PROFILES:NE_ERR", tree_name="electrons"
        )
        zts_core = params.mds_conn.get_data(
            ".YAG_NEW.RESULTS.PROFILES:Z_SORTED", tree_name="electrons"
        )
        mts_core = len(zts_core)
        zts_edge = params.mds_conn.get_data(r"\fiber_z")
        mts_edge = len(zts_edge)
        try:
            nets_edge = params.mds_conn.get_data(r"\ts_ne")
            nets_edge_err = params.mds_conn.get_data(r"\ts_ne_err")
        except mdsExceptions.mdsException as err:
            nets_edge = np.zeros((len(nets_core[:, 1]), mts_edge))
            nets_edge_err = nets_edge + 1e20
        mts = mts_core + mts_edge
        rts = params.mds_conn.get(".YAG.RESULTS.PARAM:R") + np.zeros((1, mts))
        rtci = params.mds_conn.get_data(".tci.results:rad")
        nts = len(nets_core_t)
        zts = np.zeros((1, mts))
        zts[:, :mts_core] = zts_core
        zts[:, mts_core:] = zts_edge
        nets = np.zeros((nts, mts))
        nets_err = np.zeros((nts, mts))
        nets[:, :mts_core] = (nets_core * core_mult).T
        nets_err[:, :mts_core] = (nets_core_err * core_mult).T
        nets[:, mts_core:] = (nets_edge * edge_mult).T
        nets_err[:, mts_core:] = (nets_edge_err * edge_mult).T
        valid_indices = np.where((nets_core_t >= t1) & (nets_core_t <= t2))
        if len(valid_indices) == 0:
            return t, z, n_e, n_e_sig
        nets_core_t = nets_core_t[valid_indices]
        nets = nets[valid_indices]
        nets_err = nets_err[valid_indices]
        psits = ThomsonDensityMeasure.efit_rz2psi(rts, zts, nets_core_t)
        mtci = 101
        ztci = -0.4 + 0.8 * np.arange(0, mtci) / (mtci - 1)
        rtci = rtci[nlnum] + np.zeros((1, mtci))
        psitci = ThomsonDensityMeasure.efit_rz2psi(rtci, ztci, nets_core_t)
        psia = interp1(psia_t, psia, nets_core_t)
        psi_0 = interp1(psia_t, psi_0, nets_core_t)
        nts = len(nets_core_t)
        for i in range(nts):
            psits[i, :] = (psits[i, :] - psi_0[i]) / (psia[i] - psi_0[i])
            psitci[i, :] = (psitci[i, :] - psi_0[i]) / (psia[i] - psi_0[i])
        zmapped = np.zeros((nts, 2 * mts)) + 1e32
        nemapped = zmapped.copy()
        nemapped_err = zmapped.copy()
        for i in range(nts):
            index = np.argmin(psitci[i, :]) if flag else np.argmax(psitci[i, :])
            psi_val = psitci[i, index]
            for j in range(len(mts)):
                if (flag and psits[i, j] >= psi_val) or (
                    not flag and psits[i, j] <= psi_val
                ):
                    a1 = interp1(psitci[i, :index], ztci[:index], psits[i, j])
                    a2 = interp1(psitci[i, index:], ztci[index:], psits[i, j])
                    zmapped[i, np.arange(j, j + mts + 1)] = np.arange(a1, a2)
                    nemapped[i, np.arange(j, j + mts + 1)] = nets[i, j]
                    nemapped_err[i, np.arange(j, j + mts + 1)] = nets_err[i, j]
            sorted_indices = np.argsort(zmapped[i, :])
            zmapped[i, :] = zmapped[i, sorted_indices]
            nemapped[i, :] = nemapped[i, sorted_indices]
            nemapped_err[i, :] = nemapped_err[i, sorted_indices]
        z = zmapped
        n_e = nemapped
        n_e_sig = nemapped_err
        t = nets_core_t
        return t, z, n_e, n_e_sig

    # TODO: Move to utils
    @staticmethod
    def efit_rz2psi(params: ShotDataRequestParams, r, z, t, tree="analysis"):
        r = r.flatten()
        z = z.flatten()
        psi = np.full((len(r), len(t)), np.nan)
        z = safe_cast(z, "float32")  # TODO: Ask if this change is necessary
        psirz, rgrid, zgrid, times = params.mds_conn.get_data_with_dims(
            r"\efit_geqdsk:psirz", tree_name=tree, dim_nums=[0, 1, 2]
        )
        rgrid, zgrid = np.meshgrid(rgrid, zgrid)  # , indexing='ij')

        points = np.array(
            [rgrid.flatten(), zgrid.flatten()]
        ).T  # This transposes the array to shape (n, 2)
        for i, time in enumerate(t):
            # Find the index of the closest time
            time_idx = np.argmin(np.abs(times - time))
            # Extract the corresponding Psirz slice and transpose it
            Psirz = np.transpose(psirz[time_idx, :, :])
            # Perform cubic interpolation on the Psirz slice
            values = Psirz.flatten()
            try:
                psi[:, i] = sp.interpolate.griddata(
                    points, values, (r, z), method="cubic"
                )
            except:
                params.logger.warning(
                    f"Interpolation failed for efit_rz2psi time {time}"
                )

        return psi
