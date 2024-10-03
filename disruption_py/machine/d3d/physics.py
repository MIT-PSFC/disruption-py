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


class D3DPhysicsMethods:
    """
    A class to retrieve and calculate physics-related data for DIII-D.
    """

    @staticmethod
    @physics_method(columns=["time_until_disrupt"], tokamak=Tokamak.D3D)
    def get_time_until_disrupt(params: PhysicsMethodParams):
        """
        Calculate the time until the disruption for a given shot. If the shot does
        not disrupt, return NaN.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the time until disruption. If the shot does
            not disrupt, return NaN.
        """
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
        except ValueError:
            params.logger.info(
                "[Shot %s]: Failed to get H98 signal. Returning NaNs.", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
        except ValueError:
            params.logger.info(
                "[Shot %s]: Failed to get H_alpha signal. Returning NaNs.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
        return output

    @staticmethod
    @physics_method(
        columns=["p_rad", "p_nbi", "p_ech", "radiated_fraction"],
        tokamak=Tokamak.D3D,
    )
    def get_power_parameters(params: PhysicsMethodParams):
        """
        Compute the input NBI, ECH powers, radiated power measured by the bolometer array,
        and the radiated fraction for a DIII-D shot.

        References:
        -------
        - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_power_d3d.m

        Last major update by William Wei on 8/1/2024
        """

        # Get neutral beam injected power
        try:
            p_nbi, t_nbi = params.mds_conn.get_data_with_dims(
                r"\d3d::top.nb:pinj", tree_name="d3d", astype="float64"
            )
            t_nbi /= 1e3  # [ms] -> [s]
            p_nbi *= 1e3  # [KW] -> [W]
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
                    "[Shot %s]: No NBI power data found in this shot.", params.shot_id
                )
                p_nbi = np.zeros(len(params.times))
        except mdsExceptions.MdsException:
            p_nbi = np.zeros(len(params.times))
            params.logger.info("[Shot %s]: Failed to open NBI node", params.shot_id)
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())

        # Get electron cyclotron heating (ECH) power. It's point data, so it's not
        # stored in an MDSplus tree
        try:
            p_ech, t_ech = params.mds_conn.get_data_with_dims(
                r"\top.ech.total:echpwrc", tree_name="rf"
            )
            t_ech /= 1e3  # [ms] -> [s]
            if len(t_ech) > 2:
                # Sometimes, t_ech has an extra "0" value tacked on to the end.
                # This must be removed before the interpolation.
                if t_ech[-1] == 0:
                    t_ech, p_ech = t_ech[:-1], p_ech[:-1]
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
                    "[Shot %s]: No ECH power data found in this shot. Setting to zeros",
                    params.shot_id,
                )
                p_ech = np.zeros(len(params.times))
        except mdsExceptions.MdsException:
            p_ech = np.zeros(len(params.times))
            params.logger.info(
                "[Shot %s]: Failed to open ECH node. Setting to zeros", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())

        # Get ohmic power and loop voltage
        ohmic_parameters = D3DPhysicsMethods.get_ohmic_parameters(params)
        p_ohm = ohmic_parameters["p_ohm"]

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
        except mdsExceptions.MdsException:
            params.logger.info("[Shot %s]: Failed to open bolom tree.", params.shot_id)
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
        upper_channels = [f"bol_u{i+1:02d}_v" for i in range(24)]
        lower_channels = [f"bol_l{i+1:02d}_v" for i in range(24)]
        bol_channels = upper_channels + lower_channels
        bol_signals = []
        for i in range(48):
            bol_signal = params.mds_conn.get_data(
                rf"\top.raw:{bol_channels[i]}", tree_name="bolom"
            )
            bol_signals.append(bol_signal)
        bol_time = params.mds_conn.get_dims(
            rf"\top.raw:{bol_channels[0]}", tree_name="bolom"
        )[0]
        bol_time /= 1e3  # [ms] -> [s]
        a_struct = matlab_get_bolo(
            shot_id=params.shot_id,
            bol_channels=bol_channels,
            bol_prm=bol_prm,
            bol_top=bol_signals,
            bol_time=bol_time,
            drtau=smoothing_window * 1e3,
        )
        ier = 0
        for j in range(48):
            # TODO: Ask about how many valid channels are needed for proper calculation
            if a_struct.channels[j].ier == 1:
                ier = 1
                p_rad = np.full(len(params.times), np.nan)
                break
        if ier == 0:
            b_struct = matlab_power(a_struct)
            p_rad = b_struct.pwrmix  # [W]
            p_rad = interp1(a_struct.raw_time, p_rad, params.times, "linear")

        # Remove any negative values from the power data
        # TODO: Could p_ohm be negative?
        p_rad[np.isinf(p_rad)] = np.nan
        p_rad[p_rad < 0] = 0
        p_nbi[p_nbi < 0] = 0
        p_ech[p_ech < 0] = 0

        p_input = p_ohm + p_nbi + p_ech  # [W]
        rad_fraction = p_rad / p_input
        rad_fraction[np.isinf(rad_fraction)] = np.nan

        output = {
            "p_rad": p_rad,
            "p_nbi": p_nbi,
            "p_ech": p_ech,
            "radiated_fraction": rad_fraction,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["p_ohm", "v_loop"],
        tokamak=Tokamak.D3D,
    )
    def get_ohmic_parameters(params: PhysicsMethodParams):
        """
        Compute ohmic heating power and loop voltage for a DIII-D shot

        References:
        -------
        - https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_P_ohm_d3d.m

        Last major update by William Wei on 8/1/2024
        """
        # Get edge loop voltage and smooth it a bit with a median filter
        v_loop, t_v_loop = params.mds_conn.get_data_with_dims(
            f'ptdata("vloopb", {params.shot_id})', tree_name="d3d"
        )
        t_v_loop /= 1e3  # [ms] -> [s]
        v_loop = scipy.signal.medfilt(v_loop, 11)
        v_loop = interp1(t_v_loop, v_loop, params.times, "linear")
        # Get plasma current
        ip, t_ip = params.mds_conn.get_data_with_dims(
            f"ptdata('ip', {params.shot_id})", tree_name="d3d"
        )
        t_ip /= 1e3  # [ms] -> [s]

        # Alessandro Pau (JET & AUG) has given Cristina a robust routine that
        # performs time differentiation with smoothing, while preserving causality.
        # It can be useful for differentiating numerous signals such as Ip, Vloop,
        # etc.  It is called 'GSASTD'. We will use this routine in place of Matlab's
        # 'gradient' and smoothing/filtering routines for certain signals.

        # We choose a 20-point width for gsastd. This means a 10ms window for
        # ip smoothing
        dipdt_smoothed = matlab_gsastd(
            x=t_ip,
            y=ip,
            derivative_mode=1,
            width=20,
            smooth_type=3,
            ends_type=1,
            slew_rate=0,
        )
        li, t_li = params.mds_conn.get_data_with_dims(
            r"\efit_a_eqdsk:li", tree_name="_efit_tree"
        )
        t_li /= 1e3
        # Use chisq to determine which time slices are invalid
        chisq = params.mds_conn.get_data(r"\efit_a_eqdsk:chisq", tree_name="_efit_tree")
        # Filter out invalid indices of efit reconstruction
        (invalid_indices,) = np.where(chisq > 50)
        li[invalid_indices] = np.nan

        r_0, t_r0 = params.mds_conn.get_data_with_dims(
            r"\top.results.geqdsk:rmaxis", tree_name="_efit_tree"
        )  # [m], [ms]
        t_r0 /= 1e3  # [ms] -> [s]

        li = interp1(t_li, li, params.times, "linear")
        r_0 = interp1(t_r0, r_0, params.times, "linear")
        inductance = 4.0 * np.pi * 1e-7 * r_0 * li / 2  # [H]
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
        """
        Get electron density from EFIT, then compute dn_dt and Greenwald_fraction.

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_density_parameters.m
        https://github.com/MIT-PSFC/disruption-py/issues/238
        https://github.com/MIT-PSFC/disruption-py/pull/249

        Last major update by William Wei on 8/2/2024
        """
        ne, t_ne = params.mds_conn.get_data_with_dims(
            r"\density", tree_name="_efit_tree"
        )
        # If EFIT disruption tree does not contain density data,
        # then read density from BCI subtree of D3D main tree
        # TODO: Find a shot to test this logic
        if len(~np.isnan(ne)) == 0:
            ne, t_ne = params.mds_conn.get_data_with_dims(r"\denv2", tree_name="d3d")

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
        )
        dne_dt = interp1(
            t_ne,
            dne_dt,
            params.times,
            "linear",
            bounds_error=False,
        )
        try:
            ip, t_ip = params.mds_conn.get_data_with_dims(
                f"ptdata('ip', {params.shot_id})", tree_name="_efit_tree"
            )  # [A], [ms]
            t_ip = t_ip / 1.0e3  # [ms] -> [s]
            ipsign = np.sign(np.sum(ip))
            ip = interp1(t_ip, ip * ipsign, params.times, "linear")  # positive definite
            a_minor, t_a = params.mds_conn.get_data_with_dims(
                r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
            )  # [m], [ms]
            t_a = t_a / 1.0e3  # [ms] -> [s]
            a_minor = interp1(t_a, a_minor, params.times, "linear")
            with np.errstate(divide="ignore"):
                n_g = ip / 1.0e6 / (np.pi * a_minor**2)  # [MA/m^2]
                g_f = ne / n_g * 1e-20
        except (mdsExceptions.MdsException, ValueError) as e:
            # TODO: Confirm that there is a separate exception if ptdata name doesn't exist
            params.logger.info(
                "[Shot %s]: Failed to compute Greenwald fraction.", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())

            err = "operands could not be broadcast together with shapes"
            if isinstance(ValueError, e) and err not in e.args:
                raise

            g_f = [np.nan]
        return {
            "n_e": ne,
            "greenwald_fraction": g_f,
            "dn_dt": dne_dt,
        }

    @staticmethod
    @physics_method(
        columns=["n_e_rt", "greenwald_fraction_rt", "dn_dt_rt"],
        tokamak=Tokamak.D3D,
    )
    def get_rt_density_parameters(params: PhysicsMethodParams):
        """
        Get real-time electron density from EFIT, then compute the
        real-time dn_dt and Greenwald_fraction.

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_density_parameters_RT.m
        https://github.com/MIT-PSFC/disruption-py/pull/251

        Last major update by William Wei on 8/2/2024
        """
        ne_rt, t_ne_rt = params.mds_conn.get_data_with_dims(
            f"ptdata('dssdenest', {params.shot_id})", tree_name="_efit_tree"
        )  # [10^19 m^-3]
        t_ne_rt = t_ne_rt / 1.0e3  # [ms] to [s]
        ne_rt = ne_rt * 1.0e19  # [10^19 m^-3] -> [m^-3]
        dne_dt_rt = np.gradient(ne_rt, t_ne_rt)  # [m^-3/s]
        ne_rt = interp1(t_ne_rt, ne_rt, params.times, "linear")
        dne_dt_rt = interp1(t_ne_rt, dne_dt_rt, params.times, "linear")

        # Get real time ip to calculate the Greenwald density

        try:
            ip_rt, t_ip_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('ipsip', {params.shot_id})"
            )  # [MA], [ms]
            t_ip_rt = t_ip_rt / 1.0e3  # [ms] to [s]
        except mdsExceptions.MdsException:
            ip_rt, t_ip_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('ipspr15v', {params.shot_id})"
            )  # [volts; 2 V/MA], [ms]
            t_ip_rt = t_ip_rt / 1.0e3  # [ms] to [s]
            ip_rt /= 2  # [volts] to [MA]
        ip_sign = np.sign(np.sum(ip_rt))
        ip_rt = interp1(t_ip_rt, ip_rt * ip_sign, params.times, "linear")

        # Read in EFIT minor radius and timebase.  This is also needed to calculate
        # the Greenwald density limit.  However, if the minor radius data is not
        # available, use a default fixed value of 0.59 m.  (We surveyed several
        # hundred shots to determine this default value.)  Note that the efit
        # timebase data is in a node called "atime" instead of "time" (where "time"
        # does not work).

        # For the real-time (RT) signals, read from the EFITRT1 tree
        try:
            a_minor_rt, t_a_rt = params.mds_conn.get_data_with_dims(
                r"\efit_a_eqdsk:aminor", tree_name="efitrt1"
            )  # [m], [ms]
            t_a_rt = t_a_rt / 1.0e3  # [ms] -> [s]
            a_minor_rt = interp1(t_a_rt, a_minor_rt, params.times, "linear")
        except mdsExceptions.MdsException:
            a_minor_rt = 0.59 * np.ones(len(params.times))
        try:
            with np.errstate(divide="ignore"):
                n_g_rt = ip_rt / (np.pi * a_minor_rt**2)  # [MA/m^2]
                g_f_rt = ne_rt / 1.0e20 / n_g_rt
        except ValueError as e:
            params.logger.info(
                "[Shot %s]: Failed to compute Greenwald fraction rt.", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())

            err = "operands could not be broadcast together with shapes"
            if err not in e.args:
                raise
            g_f_rt = [np.nan]

        return {"n_e_rt": ne_rt, "greenwald_fraction_rt": g_f_rt, "dn_dt_rt": dne_dt_rt}

    @staticmethod
    @physics_method(
        columns=["ip", "ip_error", "dip_dt", "dipprog_dt", "power_supply_railed"],
        tokamak=Tokamak.D3D,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """
        Retrieve plasma current parameters including measured and programmed values.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'ip' : array
                Measured plasma current values interpolated to the specified times.
            - 'ip_error' : array
                Error in plasma current, defined where feedback is active.
            - 'dip_dt' : array
                Time derivative of the measured plasma current.
            - 'dipprog_dt' : array
                Time derivative of the programmed plasma current.
            - 'power_supply_railed' : array
                Indicator of whether the power supply has railed at the specified times.
        """
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
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get measured plasma current parameters",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
                    (
                        "[Shot %s]: Polarity of Ip target is not constant. "
                        "Using value at first timestep."
                    ),
                    params.shot_id,
                )
                params.logger.debug(
                    "[Shot %s]: Polarity array %s", params.shot_id, polarity
                )
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, params.times, "linear")
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, params.times, "linear")
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get programmed plasma current parameters",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get ipimode signal. Setting to NaN.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get epsoff signal. Setting to NaN.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
            "ip_prog_rt",
            "ip_error_rt",
            "dip_dt_rt",
            "dipprog_dt_rt",
        ],
        tokamak=Tokamak.D3D,
    )
    def get_rt_ip_parameters(params: PhysicsMethodParams):
        """
        Get the real-time plasma current and programmed plasma current from EFIT,
        then compute the real-time ip_error and the derivatives of all of the above signals.

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_Ip_parameters_RT.m
        https://github.com/MIT-PSFC/disruption-py/pull/254

        Last major update by William Wei on 8/5/2024
        """
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
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get measured plasma current parameters",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
                    "[Shot %s]: Polarity of Ip target is not constant."
                    " Setting to first value in array.",
                    params.shot_id,
                )
                params.logger.debug(
                    "[Shot %s]: Polarity array: %s", params.shot_id, polarity
                )
                polarity = polarity[0]
            ip_prog_rt = ip_prog_rt * polarity
            dipprog_dt_rt = np.gradient(ip_prog_rt, t_ip_prog_rt)
            ip_prog_rt = interp1(t_ip_prog_rt, ip_prog_rt, params.times, "linear")
            dipprog_dt_rt = interp1(t_ip_prog_rt, dipprog_dt_rt, params.times, "linear")
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get programmed plasma current parameters",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
        try:
            ip_error_rt, t_ip_error_rt = params.mds_conn.get_data_with_dims(
                f"ptdata('ipeecoil', {params.shot_id})", tree_name="d3d"
            )  # [MA], [ms]
            t_ip_error_rt = t_ip_error_rt / 1.0e3  # [ms] to [s]
            ip_error_rt = ip_error_rt * 1.0e6 * 0.5  # [MA] -> [A]
            ip_error_rt = interp1(t_ip_error_rt, ip_error_rt, params.times, "linear")
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get ipeecoil signal. Setting to NaN.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
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
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get ipimode signal. Setting to NaN.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            ipimode = np.full(len(params.times), np.nan)
        (feedback_off_indices,) = np.where((ipimode != 0) & (ipimode == 3))
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
            power_supply_railed = np.zeros(len(params.times))
            (railed_indices,) = np.where(np.abs(epsoff) > 0.5)
            power_supply_railed[railed_indices] = 1
            # Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
            # PCS feedback control of Ip is not being applied.  Therefore the
            # 'ip_error' parameter is undefined for these times.
            (ps_railed_indices,) = np.where(power_supply_railed != 0)
            ip_error_rt[ps_railed_indices] = np.nan
        except mdsExceptions.MdsException:
            params.logger.info(
                (
                    "[Shot %s]: Failed to get epsoff signal. "
                    "power_supply_railed will be NaN."
                ),
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
        # 'dip_dt_RT': dip_dt_rt,
        output = {
            "ip_rt": ip_rt,
            "ip_prog_rt": ip_prog_rt,
            "ip_error_rt": ip_error_rt,
            "dip_dt_rt": dip_dt_rt,
            "dipprog_dt_rt": dipprog_dt_rt,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["zcur", "zcur_normalized"],
        tokamak=Tokamak.D3D,
    )
    def get_z_parameters(params: PhysicsMethodParams):
        """
        Get the vertical position of the plasma current centroid, then
        compute the normalized values with respect to the plasma minor radius.

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_Z_error_d3d.m
        https://github.com/MIT-PSFC/disruption-py/pull/255

        Last major update by William Wei on 9/4/2024
        """
        nominal_flattop_radius = 0.59
        # Get z_cur
        z_cur, t_z_cur = params.mds_conn.get_data_with_dims(
            f"ptdata('vpszp', {params.shot_id})", tree_name="d3d"
        )
        t_z_cur = t_z_cur / 1.0e3  # [ms] -> [s]
        z_cur = z_cur / 1.0e2  # [cm] -> [m]
        z_cur = interp1(t_z_cur, z_cur, params.times, "linear")
        # Compute z_cur_norm
        try:
            a_minor, t_a = params.mds_conn.get_data_with_dims(
                r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
            )  # [m], [ms]
            t_a = t_a / 1.0e3  # [ms] -> [s]
            chisq = params.mds_conn.get_data(
                r"\efit_a_eqdsk:chisq", tree_name="_efit_tree"
            )
            (invalid_indices,) = np.where(chisq > 50)
            a_minor[invalid_indices] = np.nan
            a_minor = interp1(t_a, a_minor, params.times, "linear")
            z_cur_norm = z_cur / a_minor
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get efit parameters", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            z_cur_norm = z_cur / nominal_flattop_radius
        return {"zcur": z_cur, "zcur_normalized": z_cur_norm}

    @staticmethod
    @physics_method(columns=["n1rms", "n1rms_normalized"], tokamak=Tokamak.D3D)
    def get_n1rms_parameters(params: PhysicsMethodParams):
        """
        Get the n1rms data, then compute n1rms_normalized = n1rms / btor

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_n1rms_d3d.m
        https://github.com/MIT-PSFC/disruption-py/pull/257

        Last major update by William Wei on 8/6/2024
        """
        # Get n1rms signal from d3d tree
        n1rms, t_n1rms = params.mds_conn.get_data_with_dims(r"\n1rms", tree_name="d3d")
        n1rms *= 1.0e-4  # Gauss -> Tesla
        t_n1rms /= 1e3  # [ms] -> [s]
        n1rms = interp1(t_n1rms, n1rms, params.times)
        # Calculate n1rms_norm
        try:
            b_tor, t_b_tor = params.mds_conn.get_data_with_dims(
                f"ptdata('bt', {params.shot_id})", tree_name="d3d"
            )
            t_b_tor /= 1e3  # [ms] -> [s]
            b_tor = interp1(t_b_tor, b_tor, params.times)  # [T]
            n1rms_norm = n1rms / np.abs(b_tor)
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to get b_tor signal to compute n1rms_normalized",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            n1rms_norm = [np.nan]
        return {"n1rms": n1rms, "n1rms_normalized": n1rms_norm}

    # TODO: Need to test and unblock recalculating peaking factors
    # By default get_peaking_factors should grab the data from MDSPlus as opposed
    # to recalculate. See DPP v4 document for details:
    # https://docs.google.com/document/d/1R7fI7mCOkMQGt8xX2nS6ZmNNkcyvPQ7NmBfRPICFaFs/edit?usp=sharing
    @staticmethod
    @physics_method(
        columns=[
            "te_peaking_cva_rt",
            "ne_peaking_cva_rt",
            "prad_peaking_cva_rt",
            "prad_peaking_xdiv_rt",
        ],
        tokamak=Tokamak.D3D,
    )
    def get_peaking_factors(params: PhysicsMethodParams):
        """
        This function calculates peaking factors for the shot number
        given by the user corresponding to the times in the given timebase.
        Electron temperature (Te_PF) and density (ne_PF) profile peaking
        factors are taken from Thomson scattering measurements, and the peaking
        factors describing radiated power distributions (Rad_CVA and Rad_XDIV)
        are taken from the 2pi foil bolometer system.

        The Thomson-based peaking factors are computed by first mapping the channel
        locations to the EFIT grid (rhovn: normalized rho, psin: normalized poloidal
        flux) and then determining the core channels through a threshold on rhovn.

        For the bolometer-based peaking factors, a subset of 12 chords from the lower
        fan array (fan = 'custom') are selected for the calculation. The core chords
        are determined through a threshold from the magnetic axis. The divertor chords
        preselected and consist of 5 chords from the 12-chord array.

        Returns
        -------
        te_peaking_cva_rt: np.ndarray
            Te peaking factor, core vs all channels
        ne_peaking_cva_rt: np.ndarray
            ne peaking factor, core vs all channels
        prad_peaking_cva_rt: np.ndarray
            bolometer peaking factor, core vs all-but-divertor channels
        prad_peaking_xdiv_rt: np.ndarray
            bolometer peaking factor, divertor vs all-but-core channels

        Reference
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_peaking_factors_d3d.m
        https://github.com/MIT-PSFC/disruption-py/pull/265
        https://github.com/MIT-PSFC/disruption-py/pull/328

        Last major update by William Wei on 10/01/2024
        """
        ## Thomson parameters
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

        ## Bolometer parameters
        # fan to use for P_rad peaking factors (either 'lower', 'upper', or 'custom')
        bolometer_fan = "custom"
        # array of bolometer fan channel numbers covering divertor
        # (upper fan: 0->23, lower fan: 24:47)
        div_channels = np.arange(26, 31)
        # time window for filtering raw bolometer signal in [ms]
        smoothing_window = 40
        p_rad_core_def = (
            0.06  # percentage of DIII-D veritcal extent defining the core margin
        )
        # 'brightness'; % either 'brightness' or 'power' ('z')
        p_rad_metric = "brightness"

        ## Additional parameters (not in MATLAB script)
        # Ts options
        ts_options = ["combined", "core", "tangential"]
        # vertical range of the DIII-D cross section in meters (for p_rad)
        vert_range = 3.0

        ne_pf = [np.nan]
        te_pf = [np.nan]
        rad_cva = [np.nan]
        rad_xdiv = [np.nan]
        # Get precomputed rad_cva & rad_xdiv data stored in ptdata tree
        calculate_prad_pf = False
        try:
            rad_cva, t_rad_cva = params.mds_conn.get_data_with_dims(
                f"ptdata('dpsrrdcva', {params.shot_id})", tree_name="d3d"
            )  # [], [ms]
            t_rad_cva /= 1e3  # [ms] -> [s]
            rad_cva = interp1(t_rad_cva, rad_cva, params.times)
            rad_xdiv, t_rad_xdiv = params.mds_conn.get_data_with_dims(
                f"ptdata('dpsrrdxdiv', {params.shot_id})", tree_name="d3d"
            )  # [], [ms]
            t_rad_xdiv /= 1e3  # [ms] -> [s]
            rad_xdiv = interp1(t_rad_xdiv, rad_xdiv, params.times)
        except mdsExceptions.MdsException:
            calculate_prad_pf = True
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            params.logger.info(
                (
                    "[Shot %s]: Failed to get rad_cva and rad_xdiv from MDSplus."
                    " Calculating using raw bolometer data."
                ),
                params.shot_id,
            )

        # Get raw Thomson data
        try:
            ts = D3DPhysicsMethods._get_ne_te(params, data_source=ts_data_type)
            for option in ts_options:
                if option in ts:
                    ts = ts[option]
                    break
            efit_dict = D3DPhysicsMethods._get_efit_dict(params)
        except (NotImplementedError, CalculationError, mdsExceptions.MdsException):
            ts = {}
            params.logger.info("[Shot %s]: Failed to get TS data", params.shot_id)
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
        if ts:
            ts["psin"], ts["rhovn"] = D3DPhysicsMethods.efit_rz_interp(ts, efit_dict)
            ts["rhovn"] = ts["rhovn"].T
            ts["psin"] = ts["psin"].T

        # Get P_rad data
        p_rad = {}
        if calculate_prad_pf:
            try:
                p_rad = D3DPhysicsMethods._get_p_rad(
                    params, fan=bolometer_fan, smoothing_window=smoothing_window
                )
            except mdsExceptions.MdsException:
                params.logger.info(
                    "[Shot %s]: Failed to get bolometer data", params.shot_id
                )
                params.logger.debug(
                    "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                )

        # Calculate te_pf & ne_pf
        if ts_radius in ts:
            # Drop data outside of valid range
            invalid_indices = np.where(
                (ts[ts_radius] < ts_radial_range[0])
                | (ts[ts_radius] > ts_radial_range[1])
            )
            ts["te"][invalid_indices] = np.nan
            ts["ne"][invalid_indices] = np.nan
            ts["te"][np.isnan(ts[ts_radius])] = np.nan
            ts["ne"][np.isnan(ts[ts_radius])] = np.nan

            # Interpolate onto uniform radial base if needed
            if ts_equispaced:
                for i in range(len(ts["time"])):
                    (no_nans,) = np.where(
                        ~np.isnan(ts["te"][:, i]) & ~np.isnan(ts["ne"][:, i])
                    )
                    if len(no_nans) <= 1:
                        continue
                    radii = ts[ts_radius][no_nans, i]
                    if len(radii) <= 2:
                        continue
                    rad_coord_interp = np.linspace(min(radii), max(radii), len(radii))
                    # MATLAB used interp1(kind='pchip') which isn't available in disruption-py
                    ts["te"][no_nans, i] = interp1(
                        radii,
                        ts["te"][no_nans, i],
                        rad_coord_interp,
                        "linear",
                    )
                    ts["ne"][no_nans, i] = interp1(
                        radii,
                        ts["ne"][no_nans, i],
                        rad_coord_interp,
                        "linear",
                    )
                    ts[ts_radius][no_nans, i] = rad_coord_interp

            # Find core bin for Thomson and calculate Te, ne peaking factors
            core_mask = ts[ts_radius] < ts_core_margin
            te_core = ts["te"].copy()
            te_core[~core_mask] = np.nan
            ne_core = ts["ne"].copy()
            ne_core[~core_mask] = np.nan
            te_pf = np.full(len(ts["time"]), np.nan)
            ne_pf = np.full(len(ts["time"]), np.nan)
            # pylint: disable-next=consider-using-enumerate
            for i in range(len(te_pf)):
                if (
                    ~np.isnan(te_core[:, i]).all()
                    and ~np.isnan(ts["te"][:, i]).all()
                    and np.nanmean(ts["te"][:, i]) != 0
                ):
                    te_pf[i] = np.nanmean(te_core[:, i]) / np.nanmean(ts["te"][:, i])
                if (
                    ~np.isnan(ne_core[:, i]).all()
                    and ~np.isnan(ts["ne"][:, i]).all()
                    and np.nanmean(ts["ne"][:, i]) != 0
                ):
                    ne_pf[i] = np.nanmean(ne_core[:, i]) / np.nanmean(ts["ne"][:, i])
            te_pf = interp1(ts["time"], te_pf, params.times)
            ne_pf = interp1(ts["time"], ne_pf, params.times)

        # Calculate prad_cva, prad_xdiv
        if calculate_prad_pf and p_rad:
            # Interpolate zmaxis and channel intersects x onto the bolometer timebase
            z_m_axis = interp1(efit_dict["time"], efit_dict["zmaxis"], p_rad["t"])
            z_m_axis = np.repeat(z_m_axis[:, np.newaxis], p_rad["x"].shape[1], axis=1)
            # NOTE: MATLAB uses extrapolation in p_rad["xinterp"] computation.
            p_rad["xinterp"] = interp1(p_rad["xtime"], p_rad["x"], p_rad["t"], axis=0)
            # Determine the bolometer channels falling in the 'core' bin
            core_indices = (
                p_rad["xinterp"] < z_m_axis + p_rad_core_def * vert_range
            ) & (p_rad["xinterp"] > z_m_axis - p_rad_core_def * vert_range)
            # Designate the divertor bin and find all 'other' channels not in that bin
            div_indices = np.full(len(p_rad["ch_avail"]), False)
            for div_channel in div_channels:
                div_indices[p_rad["ch_avail"].index(div_channel)] = True

            # Grab p_rad measurements for each needed set of channels
            p_rad_core = np.array(p_rad[p_rad_metric]).T
            p_rad_all_but_core = p_rad_core.copy()
            p_rad_div = p_rad_core.copy()
            p_rad_all_but_div = p_rad_core.copy()
            p_rad_core[~core_indices] = np.nan
            p_rad_all_but_core[core_indices] = np.nan
            p_rad_div[:, ~div_indices] = np.nan
            p_rad_all_but_div[:, div_indices] = np.nan

            # Calculate the peaking factors
            rad_cva = np.full(len(p_rad["t"]), np.nan)
            rad_xdiv = np.full(len(p_rad["t"]), np.nan)
            # pylint: disable-next=consider-using-enumerate
            for i in range(len(rad_cva)):
                if (
                    ~np.isnan(p_rad_core[i, :]).all()
                    and ~np.isnan(p_rad_all_but_div[i, :]).all()
                    and np.nanmean(p_rad_all_but_div[i, :]) != 0
                ):
                    # NOTE: How is this core vs all?
                    rad_cva[i] = np.nanmean(p_rad_core[i, :]) / np.nanmean(
                        p_rad_all_but_div[i, :]
                    )
                if (
                    ~np.isnan(p_rad_div[i, :]).all()
                    and ~np.isnan(p_rad_all_but_core[i, :]).all()
                    and np.nanmean(p_rad_all_but_core[i, :]) != 0
                ):
                    # NOTE: How is this div vs all?
                    rad_xdiv[i] = np.nanmean(p_rad_div[i, :]) / np.nanmean(
                        p_rad_all_but_core[i, :]
                    )
            rad_cva = interp1(p_rad["t"], rad_cva, params.times)
            rad_xdiv = interp1(p_rad["t"], rad_xdiv, params.times)

        output = {
            "te_peaking_cva_rt": te_pf,
            "ne_peaking_cva_rt": ne_pf,
            "prad_peaking_cva_rt": rad_cva,
            "prad_peaking_xdiv_rt": rad_xdiv,
        }
        return output

    @staticmethod
    def efit_rz_interp(ts, efit_dict):
        """
        Interpolate the efit data to the given timebase and project onto the
        poloidal plane.

        Parameters
        ----------
        ts: dict
            Thomson scattering data returned by D3DPhysicsMethods._get_ne_te(...)
        efit_dict: dict
            Dictionary with the efit data. Keys are 'time', 'r', 'z', 'psin', 'rhovn'

        Returns
        -------
        psin: np.ndarray
            Array of plasma normalized flux
        rho_vn_diag: np.ndarray
            Array of normalized minor radius

        Reference
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/sorting/efit_Rz_interp.m
        https://github.com/MIT-PSFC/disruption-py/pull/265#issuecomment-2318294825

        Last major update by William Wei on 8/29/2024
        """

        t = np.tile(ts["time"], [len(ts["r"]), 1]).transpose()
        r = np.tile(ts["r"], [len(ts["time"]), 1])
        z = np.tile(ts["z"], [len(ts["time"]), 1])

        # Implement a 3D (time,radial,vertical) gridded interpolation
        # efit_dict['psin'] has the dimensions (time, z, r)
        interp = scipy.interpolate.RegularGridInterpolator(
            [efit_dict["time"], efit_dict["z"], efit_dict["r"]],
            efit_dict["psin"],
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

        # Apply interpolant to diagnostic data and return outputs as a new structure field
        psin = interp((t, z, r))

        # Get rhovn using the interpolant stored in EFIT, then save this as another field in 'data'
        rho_vn_diag_almost = interp1(
            efit_dict["time"], efit_dict["rhovn"], ts["time"], axis=0
        )
        rho_vn_diag = np.empty(psin.shape[:2])
        # Ger the implied psin grid for rhovn
        psin_interp = np.linspace(0, 1, efit_dict["rhovn"].shape[1])
        # Interpolate again to get rhovn on same psin base
        for i in range(psin.shape[0]):
            rho_vn_diag[i] = interp1(psin_interp, rho_vn_diag_almost[i, :], psin[i, :])
        return psin, rho_vn_diag

    @staticmethod
    @physics_method(columns=["z_eff"], tokamak=Tokamak.D3D)
    def get_zeff_parameters(params: PhysicsMethodParams):
        """
        Retrieve the effective charge (Z_eff) parameters for a given shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following key:
            - 'z_eff' : array
                Effective charge values interpolated to the specified times.
        """
        # Get Zeff
        zeff, t_zeff = params.mds_conn.get_data_with_dims(
            r"\d3d::top.spectroscopy.vb.zeff:zeff", tree_name="d3d"
        )
        t_zeff = t_zeff / 1.0e3  # [ms] -> [s]
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
                "[Shot %s]: No zeff data found in this shot.", params.shot_id
            )
        return {"z_eff": zeff}

    @staticmethod
    @physics_method(columns=["kappa_area"], tokamak=Tokamak.D3D)
    def get_kappa_area(params: PhysicsMethodParams):
        """
        Compute kappa_area (elongation parameter) defined as
        plasma area / (pi * aminor**2)

        Note: the EFIT-computed kappa is retrieved in D3DEfitMethods.

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_kappa_area.m
        https://github.com/MIT-PSFC/disruption-py/pull/256

        Last major update by William Wei on 8/6/2024
        """
        a_minor = params.mds_conn.get_data(
            r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
        )
        area = params.mds_conn.get_data(r"\efit_a_eqdsk:area", tree_name="_efit_tree")
        chisq = params.mds_conn.get_data(r"\efit_a_eqdsk:chisq", tree_name="_efit_tree")
        t = params.mds_conn.get_data(r"\efit_a_eqdsk:atime", tree_name="_efit_tree")
        t /= 1e3  # [ms] -> [s]
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
        """
        Get the plasma triangularity (delta), squareness, and minor radius [m] from EFIT.

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/get_shape_parameters.m
        https://github.com/MIT-PSFC/disruption-py/pull/258

        Last major update by William Wei on 8/6/2024
        """
        # Get efit_time
        efit_time = params.mds_conn.get_data(
            r"\efit_a_eqdsk:atime", tree_name="_efit_tree"
        )
        efit_time /= 1e3  # [ms] -> [s]
        # Compute triangularity
        try:
            tritop = params.mds_conn.get_data(
                r"\efit_a_eqdsk:tritop", tree_name="_efit_tree"
            )  # meters
            tribot = params.mds_conn.get_data(
                r"\efit_a_eqdsk:tribot", tree_name="_efit_tree"
            )  # meters
            delta = (tritop + tribot) / 2.0
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to obtain triangularity signals", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            delta = [np.nan]
        # Compute squareness
        try:
            sqfod = params.mds_conn.get_data(
                r"\efit_a_eqdsk:sqfod", tree_name="_efit_tree"
            )
            sqfou = params.mds_conn.get_data(
                r"\efit_a_eqdsk:sqfou", tree_name="_efit_tree"
            )
            squareness = (sqfod + sqfou) / 2.0
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to obtain squareness signals", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            squareness = [np.nan]
        # Get aminor
        try:
            aminor = params.mds_conn.get_data(
                r"\efit_a_eqdsk:aminor", tree_name="_efit_tree"
            )
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to obtain aminor signals", params.shot_id
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
            aminor = [np.nan]
        # Remove invalid indices
        try:
            chisq = params.mds_conn.get_data(
                r"\efit_a_eqdsk:chisq", tree_name="_efit_tree"
            )
            invalid_indices = np.where(chisq > 50)
            if ~np.isnan(delta[0]):
                delta[invalid_indices] = np.nan
            if ~np.isnan(squareness[0]):
                squareness[invalid_indices] = np.nan
            if ~np.isnan(aminor[0]):
                aminor[invalid_indices] = np.nan
        except mdsExceptions.MdsException:
            params.logger.info(
                "[Shot %s]: Failed to obtain chisq to remove unreliable time points.",
                params.shot_id,
            )
            params.logger.debug("[Shot %s]: %s", params.shot_id, traceback.format_exc())
        # Interpolate to the requested time basis
        if ~np.isnan(delta[0]):
            delta = interp1(efit_time, delta, params.times, "linear")
        if ~np.isnan(squareness[0]):
            squareness = interp1(efit_time, squareness, params.times, "linear")
        if ~np.isnan(aminor[0]):
            aminor = interp1(efit_time, aminor, params.times, "linear")
        return {"delta": delta, "squareness": squareness, "aminor": aminor}

    @staticmethod
    @cache_method
    def _get_ne_te(
        params: PhysicsMethodParams,
        data_source="blessed",
        ts_systems=None,
    ):
        """
        Retrieves DIII-D Thomson scattering data

        Inputs
        -------
        data_source: string
            "blessed", "unblessed", or "ptdata'
            ("blessed" by Thomson group)
        ts_systems: list
            default: ["core", "tangential"]

        Returns
        -------
        lasers: dict

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/utils/load_ne_Te.m

        NOTE: data_source="ptdata" has not been fully implemented; however, for now this
        option isn't used in any of the methods.

        Original method by Kevin Montes on March 2019
        Last major update by William Wei on 8/8/2024
        """
        if ts_systems is None:
            ts_systems = ["core", "tangential"]
        if data_source == "blessed":  # 'blessed' by Thomson group
            mds_path = r"\top.ts.blessed."
        elif data_source == "unblessed":
            mds_path = r"\top.ts.revisions.revision00."
        elif data_source == "ptdata":
            mds_path = r"\top.ts.blessed."  # Don't ask...I don't have the answer
            raise NotImplementedError("ptdata case not fully implemented yet")  # TODO
        else:
            raise CalculationError(f"Invalid data_source: {data_source}")

        # Account for pointname formatting change in 2017 (however using ptdata is unimplemented)
        # NOTE: "suffix" is only used if data_source="ptdata" which isn't implemented yet
        suffix = {"core": "cor", "tangential": "tan"}
        if params.shot_id < 172749:  # First shot on Sep 19, 2017
            suffix["tangential"] = "hor"

        lasers = {}
        for laser in ts_systems:
            lasers[laser] = {}
            sub_tree = f"{mds_path}{laser}"
            try:
                (t_sub_tree,) = params.mds_conn.get_dims(
                    f"{sub_tree}:temp", tree_name="electrons"
                )
                # lasers[laser]['time'] gets overwritten in the loop later
                lasers[laser]["time"] = t_sub_tree / 1.0e3  # [ms] -> [s]
            except mdsExceptions.MdsException:
                lasers[laser] = None
                params.logger.info(
                    "[Shot %s]: Failed to get %s time. Setting laser data to None.",
                    params.shot_id,
                    laser,
                )
                params.logger.debug(
                    "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                )
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
                except mdsExceptions.MdsException:
                    lasers[laser][node] = np.full(lasers[laser]["time"].shape, np.nan)
                    params.logger.info(
                        "[Shot %s]: Failed to get %s:%s(%s) data, Setting to all NaNs.",
                        params.shot_id,
                        laser,
                        name,
                        node,
                    )
                    params.logger.debug(
                        "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                    )
            # Place NaNs for broken channels
            lasers[laser]["te"][lasers[laser]["te"] == 0] = np.nan
            lasers[laser]["ne"][lasers[laser]["ne"] == 0] = np.nan
            lasers[laser]["time"] /= 1e3  # [ms] -> [s]

        # If both systems/lasers available, combine them and interpolate the data
        # from the tangential system onto the finer (core) timebase
        if "tangential" in lasers and lasers["tangential"] is not None:
            if "core" in lasers and lasers["core"] is not None:
                lasers["combined"] = {}
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
        return lasers

    @staticmethod
    @cache_method
    def _get_p_rad(params: PhysicsMethodParams, fan="custom", smoothing_window=50):
        """
        Retrieves DIII-D radiation data from the bolometer MDSplus tree

        Note: a_struct.channels[i].pwr does not exactly match the results from MATLAB
        due to the use of different filtering functions (lfilter & medfilt in Python).
        However the differences are close enough so that this isn't a major problem.

        Inputs:
        -------
        fan: str
            'upper', 'lower', or 'custom' (default)

        References:
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/DIII-D/sorting/load_Prad.m

        Original author: Kevin Montes. Date: March 2019

        Last major update by William Wei on 8/30/2024
        """
        if fan == "upper":
            fan_chans = np.arange(0, 24)
        elif fan == "lower":
            fan_chans = np.arange(24, 48)
        elif fan == "custom":
            # 1st choice (heavily cover divertor and core)
            fan_chans = np.array([3, 4, 5, 6, 7, 8, 9, 12, 14, 15, 16, 22]) + 23
        else:
            return False

        # Get bolometry data
        bol_prm, _ = params.mds_conn.get_data_with_dims(r"\bol_prm", tree_name="bolom")
        upper_channels = [f"bol_u{i+1:02d}_v" for i in range(24)]
        lower_channels = [f"bol_l{i+1:02d}_v" for i in range(24)]
        bol_channels = upper_channels + lower_channels
        bol_signals = []
        bol_times = (
            []
        )  # TODO: Decide whether to actually use all bol_times instead of just first one
        for i in range(48):
            bol_signal, bol_time = params.mds_conn.get_data_with_dims(
                rf"\top.raw:{bol_channels[i]}", tree_name="bolom"
            )
            bol_time /= 1e3  # [ms] -> [s]
            bol_signals.append(bol_signal)
            bol_times.append(bol_time)
        a_struct = matlab_get_bolo(
            params.shot_id,
            bol_channels,
            bol_prm,
            bol_signals,
            bol_times[0],
            smoothing_window,
        )
        b_struct = matlab_power(a_struct)
        r_major_axis, efit_time = params.mds_conn.get_data_with_dims(
            r"\top.results.geqdsk:rmaxis", tree_name="_efit_tree"
        )
        efit_time /= 1e3  # [ms] -> [s]
        output = {
            "ch_avail": [],
            "z": [],
            "brightness": [],
            "power": [],
            "x": np.full((len(efit_time), len(fan_chans)), np.nan),
            "xtime": efit_time,
            "t": a_struct.raw_time,
        }
        if fan != "custom":
            for i, ichan in enumerate(fan_chans):
                if a_struct.channels[ichan].ier == 0:
                    output["ch_avail"].append(ichan)
                output["x"][:, i] = a_struct.channels[ichan].z + np.tan(
                    a_struct.channels[ichan].angle * np.pi / 180.0
                ) * (r_major_axis - a_struct.channels[ichan].r)
                b_struct.chan[ichan].chanpwr[
                    np.where(b_struct.chan[ichan].chanpwr < 0)
                ] = 0
                b_struct.chan[ichan].brightness[
                    np.where(b_struct.chan[ichan].brightness < 0)
                ] = 0
                output["z"].append(b_struct.chan[ichan].chanpwr)
                output["brightness"].append(b_struct.chan[ichan].brightness)
            output["power"] = output["z"]
        else:
            # All custom channels are in the lower array
            lower_fan_chans = np.arange(24, 48)
            j = 0
            for i, lower_fan_chan in enumerate(lower_fan_chans):
                # Why include these extra channels in output['power']?
                output["power"].append(b_struct.chan[lower_fan_chan].chanpwr)
                if lower_fan_chan in fan_chans:
                    ichan = fan_chans[j]
                    if a_struct.channels[ichan].ier == 0:
                        output["ch_avail"].append(ichan)
                    output["x"][:, j] = a_struct.channels[ichan].z + np.tan(
                        a_struct.channels[ichan].angle * np.pi / 180.0
                    ) * (r_major_axis - a_struct.channels[ichan].r)
                    b_struct.chan[ichan].chanpwr[
                        np.where(b_struct.chan[ichan].chanpwr < 0)
                    ] = 0
                    b_struct.chan[ichan].brightness[
                        np.where(b_struct.chan[ichan].brightness < 0)
                    ] = 0
                    output["z"].append(b_struct.chan[ichan].chanpwr)
                    output["brightness"].append(b_struct.chan[ichan].brightness)
                    j += 1
        return output

    # TODO: Replace all instances of efit_dict with a dataclass
    @staticmethod
    @cache_method
    def _get_efit_dict(params: PhysicsMethodParams):
        """
        Retrieve the EFIT data dictionary for a given shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'time' : array
                Time corresponding to the EFIT data in seconds.
            - 'z' : array
                Elevation coordinates of the grid from the EFIT data.
            - 'r' : array
                Radial coordinates of the grid from the EFIT data.
            - 'rhovn' : array
                Normalized radius from the EFIT data.
            - 'psirz' : array
                Poloidal flux on the rectangular grid points from the EFIT data.
            - 'zmaxis' : array
                Z of magnetic axis from the EFIT data.
            - 'ssimag' : array
                Poloidal flux at magnetic axis from the EFIT data.
            - 'ssibry' : array
                Poloidal flux at the plasma boundary from the EFIT data.
            - 'psin' : array
                Normalized poloidal flux values.
        """
        efit_dict = {}
        path = r"\top.results.geqdsk:"
        nodes = ["z", "r", "rhovn", "psirz", "zmaxis", "ssimag", "ssibry"]
        (efit_dict_time,) = params.mds_conn.get_dims(
            f"{path}psirz", tree_name="_efit_tree", dim_nums=[2]
        )
        efit_dict["time"] = efit_dict_time / 1e3  # [ms] -> [s]
        for node in nodes:
            try:
                efit_dict[node] = params.mds_conn.get_data(
                    f"{path}{node}", tree_name="_efit_tree"
                )
            except mdsExceptions.MdsException:
                efit_dict[node] = np.full(len(efit_dict["time"]), np.nan)
                params.logger.info(
                    "[Shot %s]: Failed to get %s from efit, Setting to all NaNs.",
                    params.shot_id,
                    node,
                )
                params.logger.debug(
                    "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                )

        # Normalize the poloidal flux grid (0=magnetic axis, 1=boundary)
        # [Translated from D. Eldon's OMFITeqdsk read_basic_eq_from_mds() function]
        psi_norm_f = efit_dict["ssibry"] - efit_dict["ssimag"]
        (problems,) = np.where(psi_norm_f == 0)
        # Prevent divide by 0 error by replacing 0s in the denominator
        psi_norm_f[problems] = 1
        efit_dict["psin"] = (
            efit_dict["psirz"] - efit_dict["ssimag"][:, np.newaxis, np.newaxis]
        ) / psi_norm_f[:, np.newaxis, np.newaxis]
        efit_dict["psin"][problems, :, :] = 0
        return efit_dict
