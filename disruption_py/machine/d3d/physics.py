#!/usr/bin/env python3

import traceback

import numpy as np
import scipy
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import get_bolo, gsastd, interp1, power
from disruption_py.machine.tokamak import Tokamak


class D3DPhysicsMethods:

    @staticmethod
    @physics_method(columns=["time_until_disrupt"], tokamak=Tokamak.D3D)
    def _get_time_until_disrupt(params: PhysicsMethodParams):
        if params.disrupted:
            return {"time_until_disrupt": params.disruption_time - params.times}
        return {"time_until_disrupt": [np.nan]}

    @staticmethod
    @physics_method(columns=["h98"], tokamak=Tokamak.D3D)
    def get_h98(params: PhysicsMethodParams):
        """
        Get the H98y2 energy confinement time parameter

        Reference
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_H98_d3d.m

        Last major update by William Wei on 7/31/2024
        """
        output = {
            "h98": [np.nan],
        }
        try:
            h_98, t_h_98 = params.mds_conn.get_data_with_dims(
                r"\H_THH98Y2", tree_name="transport"
            )
            t_h_98 /= 1e3  # [ms] -> [s]
            h_98 = interp1(t_h_98, h_98, params.times, "linear")
            output["h98"] = h_98
        except ValueError as e:
            params.logger.info(
                f"[Shot {params.shot_id}]: Failed to get H98 signal. Returning NaNs."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        return output

    @staticmethod
    @physics_method(columns=["h_alpha"], tokamak=Tokamak.D3D)
    def get_h_alpha(params: PhysicsMethodParams):
        """
        Get the H_alpha line emission intensity.

        Reference
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_H98_d3d.m

        Last major update by William Wei on 7/31/2024
        """
        output = {
            "h_alpha": [np.nan],
        }
        try:
            h_alpha, t_h_alpha = params.mds_conn.get_data_with_dims(
                r"\fs04", tree_name="d3d"
            )
            t_h_alpha /= 1e3  # [ms] -> [s]
            h_alpha = interp1(t_h_alpha, h_alpha, params.times, "linear")
            output["h_alpha"] = h_alpha
        except ValueError as e:
            params.logger.info(
                f"[Shot {params.shot_id}]: Failed to get H_alpha signal. Returning NaNs."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        return output

    @staticmethod
    @physics_method(
        columns=["p_rad", "p_nbi", "p_ech", "p_ohm", "radiated_fraction", "v_loop"],
        tokamak=Tokamak.D3D,
    )
    def get_power_parameters(params: PhysicsMethodParams):
        # Get neutral beam injected power
        try:
            p_nbi, t_nbi = params.mds_conn.get_data_with_dims(
                r"\d3d::top.nb:pinj", tree_name="d3d", astype="float64"
            )
            p_nbi *= 1.0e3  # [KW] -> [W]
            if len(t_nbi) > 2:
                p_nbi = interp1(
                    t_nbi,
                    p_nbi,
                    params.times,
                    "linear",
                    bounds_error=False,
                    fill_value=0.0,
                )
            else:
                params.logger.info(
                    f"[Shot {params.shot_id}]:No NBI power data found in this shot."
                )
                p_nbi = np.zeros(len(params.times))
        except mdsExceptions.MdsException as e:
            p_nbi = np.zeros(len(params.times))
            params.logger.info(f"[Shot {params.shot_id}]:Failed to open NBI node")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Get electron cycholotrn heating (ECH) power. It's poitn data, so it's not
        #  stored in an MDSplus tree
        try:
            p_ech, t_ech = params.mds_conn.get_data_with_dims(
                r"\top.ech.total:echpwrc", tree_name="rf"
            )
            if len(t_ech) > 2:
                p_ech = interp1(
                    t_ech,
                    p_ech,
                    params.times,
                    "linear",
                    bounds_error=False,
                    fill_value=0.0,
                )
            else:
                params.logger.info(
                    f"[Shot {params.shot_id}]:No ECH power data found in this "
                    + "shot. Setting to zeros"
                )
                p_ech = np.zeros(len(params.times))
        except mdsExceptions.MdsException as e:
            p_ech = np.zeros(len(params.times))
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to open ECH node. Setting to zeros"
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Get ohmic power and loop voltage
        p_ohm, v_loop = D3DPhysicsMethods.get_ohmic_parameters(params)
        # Radiated power
        # We had planned to use the standard signal r'\bolom::prad_tot' for this
        # parameter.  However, the processing involved in calculating \prad_tot
        # from the arrays of bolometry channels involves non-causal filtering with
        # a 50 ms window.  This is not acceptable for our purposes.  Tony Leonard
        # provided us with the two IDL routines that are used to do the automatic
        # processing that generates the \prad_tot signal in the tree (getbolo.pro
        # and powers.pro).  I converted them into Matlab routines, and modified the
        # analysis so that the smoothing is causal, and uses a shorter window.
        smoothing_window = 0.010  # [s]
        try:
            bol_prm, _ = params.mds_conn.get_data_with_dims(
                r"\bol_prm", tree_name="bolom"
            )
        except mdsExceptions.MdsException as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to open bolom tree.")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        lower_channels = [f"bol_u{i+1:02d}_v" for i in range(24)]
        upper_channels = [f"bol_l{i+1:02d}_v" for i in range(24)]
        bol_channels = lower_channels + upper_channels
        bol_signals = []
        bol_times = []
        for i in range(48):
            bol_signal, bol_time = params.mds_conn.get_data_with_dims(
                rf"\top.raw:{bol_channels[i]}", tree_name="bolom"
            )
            bol_signals.append(bol_signal)
            bol_times.append(bol_time)
        a_struct = get_bolo(
            params.shot_id, bol_channels, bol_prm, bol_signals, bol_times
        )
        ier = 0
        for j in range(48):
            # TODO: Ask about how many valid channels are needed for proper calculation
            if a_struct.channels[j].ier == 1:
                ier = 1
                p_rad = np.full(len(params.times), np.nan)
                break
        if ier == 0:
            b_struct = power(a_struct)
            p_rad = b_struct.pwrmix  # [W]
            p_rad = interp1(a_struct.time, p_rad, params.times, "linear")

        # Remove any negative values from the power data
        p_rad[np.isinf(p_rad)] = np.nan
        p_rad[p_rad < 0] = 0
        p_nbi[p_nbi < 0] = 0
        p_ech[p_ech < 0] = 0

        p_input = p_rad + p_nbi + p_ech  # [W]
        rad_fraction = p_rad / p_input
        rad_fraction[np.isinf(rad_fraction)] = np.nan

        # Computer P_sol, defined as P_in - P_rad
        p_sol = p_input - p_rad
        output = {
            "p_rad": p_rad,
            "p_nbi": p_nbi,
            "p_ech": p_ech,
            "p_ohm": p_ohm,
            "radiated_fraction": rad_fraction,
            "v_loop": v_loop,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["p_rad", "p_nbi", "p_ech", "p_ohm", "radiated_fraction", "v_loop"],
        tokamak=Tokamak.D3D,
    )
    def get_ohmic_parameters(params: PhysicsMethodParams):
        output = {
            "p_ohm": [np.nan],
            "v_loop": [np.nan],
        }
        # Get edge loop voltage and smooth it a bit with a median filter
        try:
            v_loop, t_v_loop = params.mds_conn.get_data_with_dims(
                f'ptdata("vloopb", {params.shot_id})', tree_name="d3d"
            )
            v_loop = scipy.signal.medfilt(v_loop, 11)
            v_loop = interp1(t_v_loop, v_loop, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to open VLOOPB node. Setting to NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Get plasma current
        try:
            ip, t_ip = params.mds_conn.get_data_with_dims(
                f"ptdata('ip', {params.shot_id})", tree_name="d3d"
            )
            t_ip = t_ip / 1.0e3  # [ms] -> [s]
            # We choose a 20-point width for gsastd. This means a 10ms window for
            #  ip smoothing
            dipdt_smoothed = gsastd(t_ip, ip, 1, 20, 3, 1, 0)
            li, t_li = params.mds_conn.get_data_with_dims(
                r"\efit_a_eqdsk:li", tree_name="_efit_tree"
            )
            chisq = params.mds_conn.get_data(r"\efit_a_eqdsk:chisq")
            # Filter out invalid indices of efit reconstruction
            invalid_indices = None  # TODO: Finish
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Unable to get plasma current data. p_ohm set to NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            return output
        # [m] For simplicity, use fixed r_0 = 1.67 for DIII-D major radius
        r_0 = 1.67
        inductance = 4.0 * np.pi * r_0 * li / 2  # [H]
        inductance = interp1(t_li, inductance, params.times, "linear")
        ip = interp1(t_ip, ip, params.times, "linear")
        dipdt_smoothed = interp1(t_ip, dipdt_smoothed, params.times, "linear")
        v_inductive = inductance * dipdt_smoothed  # [V]
        v_resistive = v_loop - v_inductive  # [V]
        p_ohm = ip * v_resistive  # [W]
        output = {"p_ohm": p_ohm, "v_loop": v_loop}
        return output

    @staticmethod
    @physics_method(
        columns=["n_e", "greenwald_fraction", "dn_dt"],
        tokamak=Tokamak.D3D,
    )
    def get_density_parameters(params: PhysicsMethodParams):
        ne = [np.nan]
        g_f = [np.nan]
        dne_dt = [np.nan]
        try:
            ne, t_ne = params.mds_conn.get_data_with_dims(
                r"\density", tree_name="_efit_tree"
            )
            ne = ne * 1.0e6  # [cm^3] -> [m^3]
            t_ne = t_ne / 1.0e3  # [ms] -> [s]
            dne_dt = np.gradient(ne, t_ne)
            # NOTE: t_ne has higher resolution than efit_time so t_ne[0] < efit_time[0]
            # because of rounding, meaning we need to allow extrapolation
            ne = interp1(
                t_ne,
                ne,
                params.times,
                "linear",
                bounds_error=False,
                fill_value="extrapolate",
            )
            dne_dt = interp1(
                t_ne,
                dne_dt,
                params.times,
                "linear",
                bounds_error=False,
                fill_value="extrapolate",
            )
            # TODO: CHECK TREE_NAME
            ip, t_ip = params.mds_conn.get_data_with_dims(
                f"ptdata('ip', {params.shot_id})", tree_name="_efit_tree"
            )  # [A], [ms]
            t_ip = t_ip / 1.0e3  # [ms] -> [s]
            ipsign = np.sign(np.sum(ip))
            ip = interp1(t_ip, ip * ipsign, params.times, "linear")
            a_minor, t_a = params.mds_conn.get_data_with_dims(
                r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
            )  # [m], [ms]
            t_a = t_a / 1.0e3  # [ms] -> [s]
            a_minor = interp1(t_a, a_minor, params.times, "linear")
            with np.errstate(divide="ignore"):
                n_g = ip / 1.0e6 / (np.pi * a_minor**2)  # [MA/m^2]
                g_f = ne / 1.0e20 / n_g  # TODO: Fill in units
        except mdsExceptions.MdsException as e:
            # TODO: Confirm that there is a separate exception if ptdata name doesn't exist
            params.logger.info(f"[Shot {params.shot_id}]:Failed to get some parameter")
            params.logger.debug(f"[Shot {params.shot_id}]::{traceback.format_exc()}")
        output = {"n_e": ne, "greenwald_fraction": g_f, "dn_dt": dne_dt}
        return output

    @staticmethod
    @physics_method(
        columns=["n_e_rt", "greenwald_fraction_rt"],
        tokamak=Tokamak.D3D,
    )
    def get_rt_density_parameters(params: PhysicsMethodParams):
        ne_rt = [np.nan]
        g_f_rt = [np.nan]
        dne_dt_rt = [np.nan]
        try:
            # TODO: CHECK TREE_NAME
            ne_rt, t_ne_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('dssdenest', {params.shot_id})"
            )
            t_ne_rt = t_ne_rt / 1.0e3  # [ms] to [s]
            ne_rt = ne_rt * 1.0e19  # [10^19 m^-3] -> [m^-3]
            dne_dt_rt = np.gradient(ne_rt, t_ne_rt)  # [m^-3/s]
            ne_rt = interp1(t_ne_rt, ne_rt, params.times, "linear")
            dne_dt_rt = interp1(t_ne_rt, dne_dt_rt, params.times, "linear")
            try:
                ip_rt, t_ip_rt = params.mds_conn.get_data_with_dims(
                    f"ptdata('ipsip', {params.shot_id})"
                )  # [MA], [ms]
                t_ip_rt = t_ip_rt / 1.0e3  # [ms] to [s]
                # TODO: look at units of ip_rt (not SA)
            except Exception as e:
                ip_rt, t_ip_rt = params.mds_conn.get_data_with_dims(
                    f"ptdata('ipspr15v', {params.shot_id})"
                )  # [MA], [ms]
                t_ip_rt = t_ip_rt / 1.0e3  # [ms] to [s]
            ip_sign = np.sign(np.sum(ip_rt))
            ip = interp1(t_ip_rt, ip_rt * ip_sign, params.times, "linear")
            a_minor_rt, t_a_rt = params.mds_conn.get_data_with_dims(
                r"\efit_a_eqdsk:aminor", tree_name="efitrt1"
            )  # [m], [ms]
            t_a_rt = t_a_rt / 1.0e3  # [ms] -> [s]
            a_minor_rt = interp1(t_a_rt, a_minor_rt, params.times, "linear")
            with np.errstate(divide="ignore"):
                n_g_rt = ip / 1.0e6 / (np.pi * a_minor_rt**2)  # [MA/m^2]
                g_f_rt = ne_rt / 1.0e20 / n_g_rt  # TODO: Fill in units
        except mdsExceptions.MdsException as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to get some parameter")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # ' dne_dt_RT': dne_dt_rt
        return {"n_e_rt": ne_rt, "greenwald_fraction_rt": g_f_rt}

    @staticmethod
    @physics_method(
        columns=["ip", "ip_error", "dip_dt", "dipprog_dt", "power_supply_railed"],
        tokamak=Tokamak.D3D,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        ip = [np.nan]
        ip_prog = [np.nan]
        dip_dt = [np.nan]
        dipprog_dt = [np.nan]
        # Fill with nans instead of using a single nan because indices are used
        ip_error = np.full(len(params.times), np.nan)
        # Get measured plasma current parameters
        try:
            ip, t_ip = params.mds_conn.get_data_with_dims(
                f"ptdata('ip', {params.shot_id})", tree_name="d3d"
            )  # [A], [ms]
            t_ip = t_ip / 1.0e3  # [ms] -> [s]
            dip_dt = np.gradient(ip, t_ip)
            ip = interp1(t_ip, ip, params.times, "linear")
            dip_dt = interp1(t_ip, dip_dt, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get measured plasma current parameters"
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Get programmed plasma current parameters
        try:
            ip_prog, t_ip_prog = params.mds_conn.get_data_with_dims(
                f"ptdata('iptipp', {params.shot_id})", tree_name="d3d"
            )  # [A], [ms]
            t_ip_prog = t_ip_prog / 1.0e3  # [ms] -> [s]
            polarity = np.unique(
                params.mds_conn.get_data(
                    f"ptdata('iptdirect', {params.shot_id})", tree_name="d3d"
                )
            )
            if len(polarity) > 1:
                params.logger.info(
                    f"[Shot {params.shot_id}]:Polarity of Ip target is not constant."
                    + "Using value at first timestep."
                )
                params.logger.debug(
                    f"[Shot {params.shot_id}]: Polarity array {polarity}"
                )
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, params.times, "linear")
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get programmed plasma current parameters"
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Now get the signal pointname 'ipimode'.  This PCS signal denotes whether
        # or not PCS is actually feedback controlling the plasma current.  There
        # are times when feedback of Ip is purposely turned off, such as during
        # electron cyclotron current drive experiments.  Here is how to interpret
        # the value of 'ipimode':
        #  0: normal Ip feedback to E-coils supplies
        #  3: almost normal Ip feedback, except that abs(Ip) > 2.5 MA
        #  Anything else: not in normal Ip feedback mode.  In this case, the
        # 'ip_prog' signal is irrelevant, and therefore 'ip_error' is not defined.
        try:
            ipimode, t_ipimode = params.mds_conn.get_data_with_dims(
                f"ptdata('ipimode', {params.shot_id})", tree_name="d3d"
            )
            t_ipimode = t_ipimode / 1.0e3  # [ms] -> [s]
            ipimode = interp1(t_ipimode, ipimode, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get ipimode signal. Setting to NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            ipimode = np.full(len(params.times), np.nan)
        feedback_on_indices = np.where((ipimode == 0) | (ipimode == 3))
        ip_error[feedback_on_indices] = (
            ip[feedback_on_indices] - ip_prog[feedback_on_indices]
        )
        # Finally, get 'epsoff' to determine if/when the E-coil power supplies have railed
        # Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
        # PCS feedback control of Ip is not being applied.  Therefore the
        # 'ip_error' parameter is undefined for these times.
        try:
            epsoff, t_epsoff = params.mds_conn.get_data_with_dims(
                f"ptdata('epsoff', {params.shot_id})", tree_name="d3d"
            )
            t_epsoff = t_epsoff / 1.0e3  # [ms] -> [s]
            # Avoid problem with simultaneity of epsoff being triggered exactly
            # on the last time sample
            t_epsoff += 0.001
            epsoff = interp1(t_epsoff, epsoff, params.times, "linear")
            railed_indices = np.where(np.abs(epsoff) > 0.5)
            power_supply_railed = np.zeros(len(params.times))
            power_supply_railed[railed_indices] = 1
            ip_error[railed_indices] = np.nan
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get epsoff signal. Setting to NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            power_supply_railed = [np.nan]
        # 'ip_prog': ip_prog,
        output = {
            "ip": ip,
            "ip_error": ip_error,
            "dip_dt": dip_dt,
            "dipprog_dt": dipprog_dt,
            "power_supply_railed": power_supply_railed,
        }
        return output

    @staticmethod
    @physics_method(
        columns=[
            "ip_rt",
            "ip_error_rt",
            "dipprog_dt_rt",
            "dipprog_dt",
            "power_supply_railed",
        ],
        tokamak=Tokamak.D3D,
    )
    def get_rt_ip_parameters(params: PhysicsMethodParams):
        params.mds_conn.open_tree("d3d")
        ip_rt = [np.nan]
        ip_prog_rt = [np.nan]
        ip_error_rt = [np.nan]
        dip_dt_rt = [np.nan]
        dipprog_dt_rt = [np.nan]
        # Get measured plasma current parameters
        # TODO: Why open d3d and not the rt efit tree?
        try:
            ip_rt, t_ip_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('ipsip', {params.shot_id})", tree_name="d3d"
            )  # [MA], [ms]
            t_ip_rt = t_ip_rt / 1.0e3  # [ms] -> [s]
            ip_rt = ip_rt * 1.0e6  # [MA] -> [A]
            dip_dt_rt = np.gradient(ip_rt, t_ip_rt)
            ip_rt = interp1(t_ip_rt, ip_rt, params.times, "linear")
            dip_dt_rt = interp1(t_ip_rt, dip_dt_rt, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get measured plasma current parameters"
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Get programmed plasma current parameters
        try:
            ip_prog_rt, t_ip_prog_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('ipsiptargt', {params.shot_id})", tree_name="d3d"
            )  # [MA], [ms]
            t_ip_prog_rt = t_ip_prog_rt / 1.0e3  # [ms] -> [s]
            ip_prog_rt = ip_prog_rt * 1.0e6 * 0.5  # [MA] -> [A]
            polarity = np.unique(
                params.mds_conn.get_data(
                    f"ptdata('iptdirect', {params.shot_id})", tree_name="d3d"
                )
            )
            if len(polarity) > 1:
                params.logger.info(
                    f"[Shot {params.shot_id}]:Polarity of Ip target is not constant."
                    + f" Setting to first value in array."
                )
                params.logger.debug(
                    f"[Shot {params.shot_id}]: Polarity array: {polarity}"
                )
                polarity = polarity[0]
            ip_prog_rt = ip_prog_rt * polarity
            dipprog_dt_rt = np.gradient(ip_prog_rt, t_ip_prog_rt)
            ip_prog_rt = interp1(t_ip_prog_rt, ip_prog_rt, params.times, "linear")
            dipprog_dt_rt = interp1(t_ip_prog_rt, dipprog_dt_rt, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get programmed plasma current parameters"
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        try:
            ip_error_rt, t_ip_error_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('ipeecoil', {params.shot_id})", tree_name="d3d"
            )  # [MA], [ms]
            t_ip_error_rt = t_ip_error_rt / 1.0e3  # [ms] to [s]
            ip_error_rt = ip_error_rt * 1.0e6 * 0.5  # [MA] -> [A]
            ip_error_rt = interp1(t_ip_error_rt, ip_error_rt, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get ipeecoil signal. Setting to NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Now get the signal pointname 'ipimode'.  This PCS signal denotes whether
        # or not PCS is actually feedback controlling the plasma current.  There
        # are times when feedback of Ip is purposely turned off, such as during
        # electron cyclotron current drive experiments.  Here is how to interpret
        # the value of 'ipimode':
        #  0: normal Ip feedback to E-coils supplies
        #  3: almost normal Ip feedback, except that abs(Ip) > 2.5 MA
        #  Anything else: not in normal Ip feedback mode.  In this case, the
        # 'ip_prog' signal is irrelevant, and therefore 'ip_error' is not defined.
        try:
            ipimode, t_ipimode = params.mds_conn.get_data_with_dims(
                f"ptdata('ipimode', {params.shot_id})", tree_name="d3d"
            )
            t_ipimode = t_ipimode / 1.0e3  # [ms] -> [s]
            ipimode = interp1(t_ipimode, ipimode, params.times, "linear")
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get ipimode signal. Setting to NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            ipimode = np.full(len(params.times), np.nan)
        feedback_off_indices = np.where((ipimode != 0) & (ipimode == 3))
        ip_error_rt[feedback_off_indices] = np.nan
        # Finally, get 'epsoff' to determine if/when the E-coil power supplies have railed
        # Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
        # PCS feedback control of Ip is not being applied.  Therefore the
        # 'ip_error' parameter is undefined for these times.
        try:
            epsoff, t_epsoff = params.mds_conn.get_data_with_dims(
                f"ptdata('epsoff', {params.shot_id})", tree_name="d3d"
            )
            t_epsoff = t_epsoff / 1.0e3  # [ms] -> [s]
            # Avoid problem with simultaneity of epsoff being triggered exactly on
            # the last time sample
            t_epsoff += 0.001
            epsoff = interp1(t_epsoff, epsoff, params.times, "linear")
            railed_indices = np.where(np.abs(epsoff) > 0.5)
            power_supply_railed = np.zeros(len(params.times))
            power_supply_railed[railed_indices] = 1
            ip_error_rt[railed_indices] = np.nan
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get epsoff signal. "
                + "power_supply_railed will be NaN."
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            power_supply_railed = [np.nan]
        # 'dip_dt_RT': dip_dt_rt,
        output = {
            "ip_rt": ip_rt,
            "ip_error_rt": ip_error_rt,
            "dipprog_dt_rt": dipprog_dt_rt,
            "power_supply_railed": power_supply_railed,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["zcur", "zcur_normalized", "z_prog", "z_error", "z_error_normalized"],
        tokamak=Tokamak.D3D,
    )
    def get_z_parameters(params: PhysicsMethodParams):
        """
        On DIII-D the plasma control system uses isoflux
        control to control the plasma shape and position.  It does
        NOT use zcur control.  Therefore, the PCS does not have a
        programmed vertical position.  This routine will now
        always return an arrays of NaN for z_prog, z_error, and
        z_error_norm.
        """
        NOMINAL_FLATTOP_RADIUS = 0.59
        z_cur = [np.nan]
        z_cur_norm = [np.nan]
        z_prog = [np.nan]
        z_error = [np.nan]
        z_error_norm = [np.nan]
        try:
            z_cur, t_z_cur = params.mds_conn.get_data_with_dims(
                f"ptdata('vpszp', {params.shot_id})", tree_name="d3d"
            )
            t_z_cur = t_z_cur / 1.0e3  # [ms] -> [s]
            z_cur = z_cur / 1.0e2  # [cm] -> [m]
            z_cur = interp1(t_z_cur, z_cur, params.times, "linear")
            try:
                a_minor, t_a = params.mds_conn.get_data_with_dims(
                    r"\efit_a_eqdsk:aminor", tree_name="d3d"
                )  # [m], [ms]
                t_a = t_a / 1.0e3  # [ms] -> [s]
                chisq = params.mds_conn.get_data(r"\efit_a_eqdsk:chisq")
                invalid_indices = np.where(chisq > 50)
                a_minor[invalid_indices] = np.nan
                a_minor = interp1(t_a, a_minor, params.times, "linear")
                z_cur_norm = z_cur / a_minor
            except mdsExceptions.MdsException as e:
                params.logger.info(
                    f"[Shot {params.shot_id}]:Failed to get efit parameters"
                )
                params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
                z_cur_norm = z_cur / NOMINAL_FLATTOP_RADIUS
        except mdsExceptions.MdsException as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to get vpszp signal")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        output = {
            "zcur": z_cur,
            "zcur_normalized": z_cur_norm,
            "z_prog": z_prog,
            "z_error": z_error,
            "z_error_normalized": z_error_norm,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["n_equal_1_normalized", "n_equal_1_mode"],
        tokamak=Tokamak.D3D,
    )
    def get_n1_bradial_parameters(params: PhysicsMethodParams):
        output = {
            "n_equal_1_normalized": [np.nan],
            "n_equal_1_mode": [np.nan],
        }
        # The following shots are missing bradial calculations in MDSplus and
        # must be loaded from a separate datafile
        if params.shot_id >= 176030 and params.shot_id <= 176912:
            raise NotImplementedError
            # TODO: Move to a folder like "/fusion/projects/disruption_warning/data"
            filename = "/fusion/projects/disruption_warning/matlab_programs/recalc.nc"
            # pylint: disable=undefined-variable
            ncid = nc.Dataset(filename, "r")
            brad = ncid.variables["dusbradial_calculated"][:]
            t_n1 = ncid.variables["times"][:] * 1.0e-3  # [ms] -> [s]
            shots = ncid.variables["shots"][:]
            shot_indices = np.where(shots == params.shot_id)
            if len(shot_indices) == 1:
                dusbradial = brad[shot_indices, :] * 1.0e-4  # [T]
            else:
                params.logger.info(
                    f"Shot {params.shot_id} not found in {filename}.  Returning NaN."
                )
                dusbradial = np.full(len(params.times), np.nan)
            ncid.close()
        # Check ONFR than DUD(legacy)
        else:
            try:
                # TODO: TREE NAME?
                dusbradial, t_n1 = params.mds_conn.get_data_with_dims(
                    f"ptdata('onsbradial', {params.shot_id})",
                    tree_name="d3d",
                )
                dusbradial = interp1(t_n1, dusbradial, params.times)
                dusbradial *= 1.0e-4  # [T]
            except mdsExceptions.MdsException as e:
                params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
                try:
                    dusbradial, t_n1 = params.mds_conn.get_data_with_dims(
                        f"ptdata('dusbradial', {params.shot_id})",
                        tree_name="d3d",
                    )
                    dusbradial = interp1(t_n1, dusbradial, params.times)
                    dusbradial *= 1.0e-4  # [T]
                except mdsExceptions.MdsException as e:
                    params.logger.info(
                        f"[Shot {params.shot_id}]:Failed to get n1 bradial signal. Returning NaN."
                    )
                    params.logger.debug(
                        f"[Shot {params.shot_id}]:{traceback.format_exc()}"
                    )
                    return output
        n_equal_1_mode = interp1(dusbradial, t_n1, params.times)
        # Get toroidal field Btor
        b_tor, t_b_tor = params.mds_conn.get_data_with_dims(
            f"ptdata('bt', {params.shot_id})", tree_name="d3d"
        )
        b_tor = interp1(t_b_tor, b_tor, params.times)  # [T]
        n_equal_1_normalized = n_equal_1_mode / b_tor
        output = {
            "n_equal_1_normalized": n_equal_1_normalized,
            "n_equal_1_mode": n_equal_1_mode,
        }
        return output

    @staticmethod
    @physics_method(columns=["n1rms", "n1rms_normalized"], tokamak=Tokamak.D3D)
    def get_n1rms_parameters(params: PhysicsMethodParams):
        """
        Get n1rms data, then compute n1rms_normalized = n1rms / btor

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_n1rms_d3d.m
        https://github.com/MIT-PSFC/disruption-py/pull/257

        Last major update by William Wei on 8/6/2024
        """
        # Get n1rms signal from d3d tree
        try:
            n1rms, t_n1rms = params.mds_conn.get_data_with_dims(
                r"\n1rms", tree_name="d3d"
            )
            n1rms *= 1.0e-4  # Gauss -> Tesla
            t_n1rms /= 1e3  # [ms] -> [s]
            n1rms = interp1(t_n1rms, n1rms, params.times)
        except mdsExceptions.MdsException as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to get n1rms signal")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            return {"n1rms": [np.nan], "n1rms_normalized": [np.nan]}
        # Calculate n1rms_norm
        try:
            b_tor, t_b_tor = params.mds_conn.get_data_with_dims(
                f"ptdata('bt', {params.shot_id})", tree_name="d3d"
            )
            t_b_tor /= 1e3  # [ms] -> [s]
            b_tor = interp1(t_b_tor, b_tor, params.times)  # [T]
            n1rms_norm = n1rms / np.abs(b_tor)
        except mdsExceptions.MdsException as e:
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get b_tor signal to compute n1rms_normalized"
            )
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            n1rms_norm = [np.nan]
        return {"n1rms": n1rms, "n1rms_normalized": n1rms_norm}

    # TODO: Need to test and unblock recalculating peaking factors
    # By default get_peaking_factors should grab the data from MDSPlus as opposed
    # to recalculate. See DPP v4 document for details:
    # https://docs.google.com/document/d/1R7fI7mCOkMQGt8xX2nS6ZmNNkcyvPQ7NmBfRPICFaFs/edit?usp=sharing
    @staticmethod
    @physics_method(
        columns=["te_pf", "ne_pf", "rad_cva", "rad_xdiv"],
        tokamak=Tokamak.D3D,
    )
    def get_peaking_factors(params: PhysicsMethodParams):
        ts_data_type = "blessed"  # either 'blessed', 'unblessed', or 'ptdata'
        # metric to use for core/edge binning (either 'psin' or 'rhovn')
        ts_radius = "rhovn"
        # ts_radius value defining boundary of 'core' region (between 0 and 1)
        ts_core_margin = 0.3
        # All data outside this range excluded. For example, psin=0 at magnetic axis
        # and 1 at separatrix.
        ts_radial_range = (0, 1)
        # set to true to interpolate ts_channel data onto equispaced radial grid
        ts_equispaced = False
        # fan to use for P_rad peaking factors (either 'lower', 'upper', or 'custom')
        bolometer_fan = "custom"
        # array of bolometer fan channel numbers covering divertor
        # (upper fan: 1->24, lower fan: 25:48)
        div_channels = np.arange(3, 8) + 24
        # time window for filtering raw bolometer signal in [ms]
        smoothing_window = 40
        p_rad_core_def = (
            0.06  # percentage of DIII-D veritcal extent defining the core margin
        )
        # 'brightness'; % either 'brightness' or 'power' ('z')
        p_rad_metric = "brightness"
        # Ts options
        ts_options = ["combined", "core", "tangential"]
        # vertical range of the DIII-D cross section in meters
        vert_range = 3.0
        ne_pf = [np.nan]
        te_pf = [np.nan]
        rad_cva = [np.nan]
        rad_xdiv = [np.nan]
        try:
            # TODO: TREE NAME
            rad_cva, t_rad_cva = params.mds_conn.get_data_with_dims(
                f"ptdata('dpsradcva', {params.shot_id})", tree_name="d3d"
            )
            rad_cva = interp1(t_rad_cva, rad_cva, params.times)  # [T]

            rad_xdiv, t_rad_xdiv = params.mds_conn.get_data_with_dims(
                f"ptdata('dpsradxdiv', {params.shot_id})", tree_name="d3d"
            )
            rad_xdiv = interp1(t_rad_xdiv, rad_xdiv, params.times)  # [T]
        except mdsExceptions.MdsException as e:
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            params.logger.info(
                f"[Shot {params.shot_id}]:Failed to get CVA and XDIV from MDSPlus."
                + " Calculating locally, results may be inaccurate."
            )
            rad_cva = [np.nan]
            rad_xdiv = [np.nan]
        try:
            ts = D3DPhysicsMethods._get_ne_te(params)
            for option in ts_options:
                if option in ts:
                    ts = ts[option]
            efit_dict = D3DPhysicsMethods._get_efit_dict(params)
        except Exception as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to get TS data")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            ts = 0
        try:
            ts["psin"], ts["rhovn"] = D3DPhysicsMethods.efit_rz_interp(ts, efit_dict)
            ts["rhovn"] = ts["rhovn"].T
            ts["psin"] = ts["psin"].T
            params.logger.info(ts["rhovn"].shape)
        except Exception as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to interpolate TS data")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        try:
            p_rad = D3DPhysicsMethods._get_p_rad(params)
        except Exception as e:
            params.logger.info(f"[Shot {params.shot_id}]:Failed to get bolometer data")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
            p_rad = 0
        if p_rad == 0 and ts == 0:
            params.logger.info(
                f"[Shot {params.shot_id}]:Both TS and bolometer data missing for shot"
            )

        # if ts_equispaced:
        if ts != 0 and ts_radius in ts:
            # Drop data outside of valid range
            invalid_indices = np.where(
                (ts[ts_radius] < ts_radial_range[0])
                | (ts[ts_radius] > ts_radial_range[1])
            )
            ts["te"][invalid_indices] = np.nan
            ts["ne"][invalid_indices] = np.nan
            ts["te"][np.isnan(ts[ts_radius])] = np.nan
            ts["ne"][np.isnan(ts[ts_radius])] = np.nan
            if ts_equispaced:
                raise NotImplementedError(
                    "Equispaced is currently assumed to be false"
                )  # TODO
            # Find core bin for Thomson and calculate Te, ne peaking factors
            core_mask = ts[ts_radius] < ts_core_margin
            te_core = ts["te"]
            te_core[~core_mask] = np.nan
            ne_core = ts["ne"]
            ne_core[~core_mask] = np.nan
            te_pf = np.nanmean(te_core, axis=0) / np.nanmean(ts["te"], axis=0)
            ne_pf = np.nanmean(ne_core, axis=0) / np.nanmean(ts["ne"], axis=0)
            te_pf = interp1(ts["time"], te_pf, params.times)
            ne_pf = interp1(ts["time"], ne_pf, params.times)
            # Calculate Prad CVA, X-DIV Peaking Factors
            # # Interpolate zmaxis and channel intersects x onto the bolometer timebase
            z_m_axis = interp1(efit_dict["time"], efit_dict["zmaxis"], p_rad["t"])
            z_m_axis = np.repeat(z_m_axis[:, np.newaxis], p_rad["x"].shape[1], axis=1)
            p_rad["xinterp"] = interp1(p_rad["xtime"], p_rad["x"], p_rad["t"], axis=0)
            # # Determine the bolometer channels falling in the 'core' bin
            core_indices = (
                p_rad["xinterp"] < z_m_axis + p_rad_core_def * vert_range
            ) & (p_rad["xinterp"] > z_m_axis - p_rad_core_def * vert_range)
            # # Designate the divertor bin and find all 'other' channels not in that bin
            div_indices = np.searchsorted(p_rad["ch_avail"], div_channels)
            other_indices = ~div_indices
            # # Grab p_rad measurements for each needed set of channels
            p_rad_core = np.array(p_rad[p_rad_metric]).T
            p_rad_all_but_core = p_rad_core.copy()
            p_rad_div = p_rad_core.copy()
            p_rad_all_but_div = p_rad_core.copy()
            # QUESTION: Why fill with nans for core but just keep valid indices for divertor
            p_rad_core[~core_indices] = np.nan
            p_rad_all_but_core[core_indices] = np.nan
            p_rad_div = p_rad_div[:, div_indices]
            p_rad_all_but_div = p_rad_all_but_div[:, other_indices]
            # # Calculate the peaking factors
            rad_cva = np.nanmean(p_rad_core, axis=1) / np.nanmean(
                p_rad_all_but_div, axis=1
            )
            rad_xdiv = np.nanmean(p_rad_div, axis=1) / np.nanmean(
                p_rad_all_but_core, axis=1
            )
            rad_cva = interp1(p_rad["t"], rad_cva.T, params.times)
            rad_xdiv = interp1(p_rad["t"], rad_xdiv.T, params.times)
        output = {
            "te_pf": te_pf,
            "ne_pf": ne_pf,
            "rad_cva": rad_cva,
            "rad_xdiv": rad_xdiv,
        }
        return output

    @staticmethod
    # TODO: Finish implementing just in case
    def _efit_map_rz_to_rho_original(params: PhysicsMethodParams, ts_dict, efit_dict):
        slices = np.zeros(ts_dict["time"].shape)
        # If thomson starts before EFIT (often does), then use the first valid EFIT
        # slice for early Thomson data.
        early_indices = np.where(ts_dict["time"] < efit_dict["time"])
        if len(early_indices[0]) > 0:
            slices[early_indices] = 1
            first_ts = early_indices[0][-1]
        else:
            first_ts = 0
        # If Thomson ends after EFIT (also often happens), then use the last valid EFIT
        # slice for late Thomson data.
        late_indices = np.where(ts_dict["time"] >= efit_dict["time"])
        if len(late_indices[0]) > 0:
            slices[late_indices] = len(efit_dict["time"])
            last_ts = late_indices[0][0] - 1
        else:
            last_ts = len(ts_dict["time"]) - 1
        diag_slices = np.arange(first_ts, last_ts + 1, 1)
        # Acquire list of diag time slices w/in EFIT time range; Should find closest EFIT
        # for each one
        for i in diag_slices:
            slices[i] = np.argmin(np.abs(efit_dict["time"] - ts_dict["time"][i]))
        # Interpolate EFIT data onto Thomson time slices
        psin_diag_arr = np.zeros((len(efit_dict["time"]), len(ts_dict["z"])))
        for r in np.unique(ts_dict["r"]):
            dr = r - efit_dict["r"]
            # Find closet EFIT R on the left and right
            right = np.where(efit_dict["r"] > r, 1)
            left = right - 1
            if efit_dict["r"][right] == r:
                psin_slice = np.squeeze(efit_dict["psin"][:, right, :])

    @staticmethod
    def efit_rz_interp(ts, efit_dict):
        """
        Interpolate the efit data to the given timebase and project onto the
        poloidal plane.
        Parameters
        ----------
        ts: np.ndarray
            Timebase to interpolate to
        efit_dict: dict
            Dictionary with the efit data. Keys are 'time', 'r', 'z', 'psin', 'rhovn'
        Returns
        -------
        psin: np.ndarray
            Array of plasma normalized flux
        rho_vn_diag: np.ndarray
            Array of normalized minor radius
        """
        times = ts["time"] / 1.0e3
        interp = scipy.interpolate.RegularGridInterpolator(
            [efit_dict["time"], efit_dict["r"], efit_dict["z"]],
            efit_dict["psin"],
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        # T,R,Z = np.meshgrid(times, efit_dict['r'], efit_dict['z'],indexing='ij')
        T, R, Z = np.meshgrid(times, ts["r"], ts["z"], indexing="ij")
        print("EFIT rhovn shape:", efit_dict["rhovn"].shape)
        # print(np.stack((T,R,Z),axis=1).shape)
        psin = interp((T, R, Z))
        rho_vn_diag_almost = interp1(
            efit_dict["time"], efit_dict["rhovn"], times, axis=0
        )
        print("Rho_vn_diag_almost shape", rho_vn_diag_almost.shape)
        rho_vn_diag = np.empty(psin.shape[:2])
        psin_timebase = np.linspace(0, 1, efit_dict["rhovn"].shape[1])
        for i in range(psin.shape[0]):
            rho_vn_diag[i] = interp1(
                psin_timebase, rho_vn_diag_almost[i,], psin[i, :]
            ).diagonal()
        return psin, rho_vn_diag

    @staticmethod
    @physics_method(columns=["z_eff"], tokamak=Tokamak.D3D)
    def get_zeff_parameters(params: PhysicsMethodParams):
        # Get Zeff
        try:
            zeff, t_zeff = params.mds_conn.get_data_with_dims(
                r"\d3d::top.spectroscopy.vb.zeff:zeff", tree_name="d3d"
            )
            t_zeff = t_zeff / 1.0e3  # [ms] -> [s]
            # t_nbi = params.mds_conn.get(
            # r"dim_of(\d3d::top.nb:pinj)").data()/1.e3  # [ms]->[s]
            if len(t_zeff) > 2:
                zeff = interp1(
                    t_zeff,
                    zeff,
                    params.times,
                    "linear",
                    bounds_error=False,
                    fill_value=0.0,
                )
            else:
                zeff = np.zeros(len(params.times))
                params.logger.info(
                    f"[Shot {params.shot_id}]:No zeff data found in this shot."
                )
        except mdsExceptions.MdsException as e:
            zeff = np.zeros(len(params.times))
            params.logger.info(f"[Shot {params.shot_id}]:Failed to open Zeff node")
            params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        return {"z_eff": zeff}

    @staticmethod
    @physics_method(columns=["kappa_area"], tokamak=Tokamak.D3D)
    def get_kappa_area(params: PhysicsMethodParams):
        a_minor = params.mds_conn.get_data(
            r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
        )
        area = params.mds_conn.get_data(r"\efit_a_eqdsk:area", tree_name="_efit_tree")
        chisq = params.mds_conn.get_data(r"\efit_a_eqdsk:chisq", tree_name="_efit_tree")
        t = params.mds_conn.get(r"\efit_a_eqdsk:atime", tree_name="_efit_tree")
        kappa_area = area / (np.pi * a_minor**2)
        invalid_indices = np.where(chisq > 50)
        kappa_area[invalid_indices] = np.nan
        kappa_area = interp1(t, kappa_area, params.times)
        return {"kappa_area": kappa_area}

    @staticmethod
    @physics_method(
        columns=["delta", "squareness", "aminor"],
        tokamak=Tokamak.D3D,
    )
    def get_shape_parameters(params: PhysicsMethodParams):
        efit_time = (
            params.mds_conn.get_data(r"\efit_a_eqdsk:atime", tree_name="_efit_tree")
            / 1.0e3
        )  # [ms] -> [s]
        sqfod = params.mds_conn.get_data(r"\efit_a_eqdsk:sqfod", tree_name="_efit_tree")
        sqfou = params.mds_conn.get_data(r"\efit_a_eqdsk:sqfou", tree_name="_efit_tree")
        tritop = params.mds_conn.get_data(
            r"\efit_a_eqdsk:tritop", tree_name="_efit_tree"
        )  # meters
        tribot = params.mds_conn.get_data(
            r"\efit_a_eqdsk:tribot", tree_name="_efit_tree"
        )  # meters
        # plasma minor radius [m]
        aminor = params.mds_conn.get_data(
            r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
        )
        chisq = params.mds_conn.get_data(r"\efit_a_eqdsk:chisq", tree_name="_efit_tree")
        # Compute triangularity and squareness:
        delta = (tritop + tribot) / 2.0
        squareness = (sqfod + sqfou) / 2.0

        # Remove invalid indices
        invalid_indices = np.where(chisq > 50)
        delta[invalid_indices] = np.nan
        squareness[invalid_indices] = np.nan
        aminor[invalid_indices] = np.nan

        # Interpolate to desired times
        delta = interp1(
            efit_time,
            delta,
            params.times,
            "linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        squareness = interp1(
            efit_time,
            squareness,
            params.times,
            "linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        aminor = interp1(
            efit_time,
            aminor,
            params.times,
            "linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        output = {"delta": delta, "squareness": squareness, "aminor": aminor}
        return output

    @staticmethod
    @cache_method
    def _get_ne_te(
        params: PhysicsMethodParams,
        data_source="blessed",
        ts_systems=["core", "tangential"],
    ):
        if data_source == "blessed":  # 'blessed' by Thomson group
            mds_path = r"\top.ts.blessed."
        elif data_source == "unblessed":
            mds_path = r"\top.ts.revisions.revision00."
        elif data_source == "ptdata":
            mds_path = r"\top.ts.blessed."  # Don't ask...I don't have the answer
            raise NotImplementedError("ptdata case not fully implemented yet")  # TODO
        else:
            raise ValueError(f"Invalid data_source: {data_source}")
        # Account for pointname formatting change in 2017 (however using ptdata is unimplemented)
        suffix = {"core": "cor", "tangential": "tan"}
        if params.shot_id < 172749:  # First shot on Sep 19, 2017
            suffix["tangential"] = "hor"
        lasers = dict()
        for laser in ts_systems:
            lasers[laser] = dict()
            sub_tree = f"{mds_path}{laser}"
            try:
                (t_sub_tree,) = params.mds_conn.get_dims(
                    f"{sub_tree}:temp", tree_name="electrons"
                )
                lasers[laser]["time"] = t_sub_tree / 1.0e3  # [ms] -> [s]
            except mdsExceptions.MdsException as e:
                lasers[laser] = None
                params.logger.info(
                    f"[Shot {params.shot_id}]: Failed to get {laser} time. Setting laser data to None."
                )
                params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
                continue
            child_nodes = {
                "r": "r",
                "z": "z",
                "te": "temp",
                "ne": "density",
                "time": "time",
                "te_error": "temp_e",
                "ne_error": "density_e",
            }
            for node, name in child_nodes.items():
                try:
                    lasers[laser][node] = params.mds_conn.get_data(
                        f"{sub_tree}:{name}", tree_name="electrons"
                    )
                except mdsExceptions.MdsException as e:
                    lasers[laser][node] = np.full(lasers[laser]["time"].shape, np.nan)
                    params.logger.info(
                        f"[Shot {params.shot_id}]: Failed to get {laser}:{name}({node})"
                        + " data, Setting to all NaNs."
                    )
                    params.logger.debug(
                        f"[Shot {params.shot_id}]:{traceback.format_exc()}"
                    )
            # Place NaNs for broken channels
            lasers[laser]["te"][lasers[laser]["te"] == 0] = np.nan
            lasers[laser]["ne"][np.where(lasers[laser]["ne"] == 0)] = np.nan
            params.logger.debug(
                "_get_ne_te: Core bins {}".format(lasers["core"]["te"].shape)
            )
            params.logger.debug(
                "_get_ne_te: Tangential bins {}".format(
                    lasers["tangential"]["te"].shape
                )
            )
        # If both systems/lasers available, combine them and interpolate the data
        # from the tangential system onto the finer (core) timebase
        if "tangential" in lasers and lasers["tangential"] is not None:
            if "core" in lasers and lasers["core"] is not None:
                lasers["combined"] = dict()
                # Interpolate tangential data onto core timebase
                for key in lasers["tangential"]:
                    if key not in ["time", "r", "z"]:
                        lasers["tangential"][key] = interp1(
                            lasers["tangential"]["time"],
                            lasers["tangential"][key],
                            lasers["core"]["time"],
                        )
                        lasers["combined"][key] = np.concatenate(
                            (lasers["core"][key], lasers["tangential"][key])
                        )
                lasers["tangential"]["time"] = lasers["core"]["time"]
                lasers["combined"]["time"] = lasers["core"]["time"]
                lasers["combined"]["r"] = np.concatenate(
                    (lasers["core"]["r"], lasers["tangential"]["r"])
                )
                lasers["combined"]["z"] = np.concatenate(
                    (lasers["core"]["z"], lasers["tangential"]["z"])
                )
        params.logger.debug("_get_ne_te: R Bins:", len(lasers["combined"]["r"]))
        params.logger.debug("_get_ne_te: Z Bins:", len(lasers["combined"]["z"]))
        return lasers

    @staticmethod
    @cache_method
    def _get_p_rad(params: PhysicsMethodParams, fan="custom"):
        if fan == "upper":
            fan_chans = np.arange(0, 24)
        elif fan == "lower":
            fan_chans = np.arange(24, 48)
        # Default is fan="custom"
        else:
            # 1st choice (heavily cover divertor and core)
            fan_chans = np.array([3, 4, 5, 6, 7, 8, 9, 12, 14, 15, 16, 22]) + 24

        # Get bolometry data
        bol_prm, _ = params.mds_conn.get_data_with_dims(r"\bol_prm", tree_name="bolom")
        lower_channels = [f"bol_u{i+1:02d}_v" for i in range(24)]
        upper_channels = [f"bol_l{i+1:02d}_v" for i in range(24)]
        bol_channels = lower_channels + upper_channels
        bol_signals = []
        bol_times = (
            []
        )  # TODO: Decide whether to actually use all bol_times instead of just first one
        for i in range(48):
            bol_signal, bol_time = params.mds_conn.get_data_with_dims(
                rf"\top.raw:{bol_channels[i]}", tree_name="bolom"
            )
            bol_signals.append(bol_signal)
            bol_times.append(bol_time)
        a_struct = get_bolo(
            params.shot_id, bol_channels, bol_prm, bol_signals, bol_times[0]
        )
        b_struct = power(a_struct)
        r_major_axis, efit_time = params.mds_conn.get_data_with_dims(
            r"\top.results.geqdsk:rmaxis", tree_name="_efit_tree"
        )
        output = {
            "ch_avail": [],
            "z": [],
            "brightness": [],
            "power": [],
            "x": np.full((len(efit_time), len(fan_chans)), np.nan),
            "xtime": efit_time,
            "t": a_struct.raw_time,
        }
        for i in range(len(fan_chans)):
            chan = fan_chans[i]
            output["power"].append(b_struct.chan[chan].chanpwr)
            if a_struct.channels[chan].ier == 0:
                output["ch_avail"].append(chan)
            output["x"][:, i] = a_struct.channels[chan].Z + np.tan(
                a_struct.channels[chan].angle * np.pi / 180.0
            ) * (r_major_axis - a_struct.channels[chan].R)
            b_struct.chan[chan].chanpwr[np.where(b_struct.chan[chan].chanpwr < 0)] = 0
            b_struct.chan[chan].brightness[
                np.where(b_struct.chan[chan].brightness < 0)
            ] = 0
            output["z"].append(b_struct.chan[i].chanpwr)
            output["brightness"].append(b_struct.chan[i].brightness)
        return output

    # TODO: Replace all instances of efit_dict with a dataclass
    @staticmethod
    @cache_method
    def _get_efit_dict(params: PhysicsMethodParams):
        efit_dict = dict()
        path = r"\top.results.geqdsk:"
        nodes = ["z", "r", "rhovn", "psirz", "zmaxis", "ssimag", "ssibry"]
        (efit_dict_time,) = params.mds_conn.get_dims(
            f"{path}psirz", tree_name="_efit_tree", dim_nums=[2]
        )
        efit_dict["time"] = efit_dict_time / 1.0e3  # [ms] -> [s]
        for node in nodes:
            try:
                efit_dict[node] = params.mds_conn.get_data(
                    f"{path}{node}", tree_name="_efit_tree"
                )
            except mdsExceptions.MdsException as e:
                efit_dict[node] = np.full(len(efit_dict["time"]), np.nan)
                params.logger.info(
                    f"[Shot {params.shot_id}]: Failed to get {node} from efit, Setting to all NaNs."
                )
                params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        # Normalize the poloidal flux grid (0=magnetic axis, 1=boundary)
        # [Translated from D. Eldon's OMFITeqdsk read_basic_eq_from_mds() function]
        psi_norm_f = efit_dict["ssibry"] - efit_dict["ssimag"]
        problems = np.where(psi_norm_f == 0)[0]
        # Prevent divide by 0 error by replacing 0s in the denominator
        psi_norm_f[problems] = 1
        efit_dict["psin"] = (
            efit_dict["psirz"] - efit_dict["ssimag"][:, np.newaxis, np.newaxis]
        ) / psi_norm_f[:, np.newaxis, np.newaxis]
        efit_dict["psin"][problems, :, :] = 0
        return efit_dict
