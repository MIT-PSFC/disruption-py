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
    smooth,
)
from disruption_py.machine.tokamak import Tokamak


class EASTPhysicsMethods:
    """
    A class to retrieve and calculate physics-related data for DIII-D.
    """

    @staticmethod
    @physics_method(columns=["time_until_disrupt"], tokamak=Tokamak.EAST)
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

        Last major update: 11/19/24 by William Wei
        """
        if params.disrupted:
            return {"time_until_disrupt": params.disruption_time - params.times}
        return {"time_until_disrupt": [np.nan]}

    @staticmethod
    @physics_method(
        columns=[
            "ip",
            "ip_prog",
            "ip_error",
            "ip_error_normalized",
            "dip_dt",
            "dipprog_dt",
        ],
        tokamak=Tokamak.EAST,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """
        This routine calculates Ip_error = (Ip - Ip_programmed), i.e. how much
        the actual plasma current differs from the requested current.  It
        linearly interpolates both the programmed and measured plasma currents
        onto the given timebase.  The routine also calculates the time
        derivatives of the measured Ip and programmed Ip.  The time derivatives
        are useful for discriminating between rampup, flattop, and rampdown.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'ip' : array
                Measured plasma current [A].
            - 'ip_prog' : array
                Programmed (requested) plasma current [A].
            - 'ip_error' : array
                ip - ip_prog [A].
            - 'ip_error_normalized' : array
                ip_error normalized to ip_prog [dimensionless].
            - 'dip_dt' : array
                Time derivative of the measured plasma current [A/s].
            - 'dipprog_dt' : array
                Time derivative of the programmed plasma current [A/s]

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_Ip_parameters.m

        Original author: Robert Granetz, Dec 2015

        Last major update: 11/19/24 by William Wei
        """
        ip = [np.nan]
        ip_prog = [np.nan]
        ip_error = [np.nan]
        ip_error_normalized = [np.nan]
        dip_dt = [np.nan]
        dipprog_dt = [np.nan]

        # Read in the measured plasma current, Ip.  There are several
        # different measurements of Ip: IPE, IPG, IPM (all in the EAST tree), and
        # PCRL01 (in the PCS_EAST tree).  At various times in the history of EAST,
        # there are been problems with all of these measurements, such as broken
        # sensors, inverted signals, and shifted timebases.  I think the most
        # reliable one is PCRL01, which is the one used by the Plasma Control
        # System (PCS) for feedback control.  So that is the one I will use for the
        # disruption warning database.
        ip, ip_time = params.mds_conn.get_data_with_dims(
            r"\pcrl01", tree_name="pcs_east"
        )  # [A], [s]

        # For shots before year 2014, the PCRL01 timebase needs to be shifted
        # by 17.0 ms
        if params.shot_id < 44432:
            ip_time -= 0.0170

        # High-frequency noise spikes on some shots can cause a problem with the
        # time derivative and other computations.  Use a median filter to reduce
        # the problem.
        ip = scipy.signal.medfilt(ip, 5)  # Remove noise spikes with median filter

        # Subtract baseline offset
        (base_indices,) = np.where(
            ip_time <= -5.8
        )  # time before any PF supplies turn on
        if len(base_indices) > 0:
            baseline = sum(ip[base_indices]) / len(base_indices)
            ip -= baseline

        # Calculate dip_dt
        dip_dt = np.gradient(ip, ip_time)

        # Get the programmed plasma current.  For most times, the EAST plasma
        # control system (PCS) programming is signal "\lmtipref" in the pcs_east
        # tree.  However, there is an alternate control scheme, "isoflux" control,
        # which uses one of four possible signals (\ietip, \idtip, \istip, iutip).
        # Note that only one of these 4 signals is defined at any given time during
        # the shot.  The transition from normal control to isoflux control is
        # identified by several signals, one of which is called "\sytps1".  When
        # \sytps1 = 0, use \lmtipref for ip_prog.  When \sytps1 > 0, use whichever
        # isoflux or transition signal is defined.  Also, accordingly to Yuan
        # Qiping (EAST PCS expert), no shot is ever programmed to have a
        # back-transition, i.e. after a transition to isoflux control, no shot is
        # ever programmed to transition back to standard control.

        # Start with \lmtipref (standard Ip programming), which is defined for the
        # entire shot, and defines the timebase for the programmed Ip signal.
        ip_prog, ip_prog_time = params.mds_conn.get_data_with_dims(
            r"\lmtipref*1e6", tree_name="pcs_east"
        )  # [A], [s]

        # Now check to see if there is a transition to isoflux control
        sytps1, time_sytps1 = params.mds_conn.get_data_with_dims(
            r"\sytps1", tree_name="pcs_east"
        )
        (sytps1_indices,) = np.where(sytps1 > 0)

        # If PCS switches to isoflux control, then read in the 4 possible Ip target
        # signals, and interpolate each one onto the \lmtipref timebase.  There will
        # be many NaN values.

        if len(sytps1_indices) > 0 and time_sytps1[0] <= ip_prog_time[-1]:
            try:
                ietip, ietip_time = params.mds_conn.get_data_with_dims(
                    r"\ietip", tree_name="pcs_east"
                )
                ietip = interp1(ietip_time, ietip, ip_prog_time)
            except mdsExceptions.MdsException:
                ietip = np.full(len(ip_prog_time), np.nan)

            try:
                idtip, idtip_time = params.mds_conn.get_data_with_dims(
                    r"\idtip", tree_name="pcs_east"
                )
                idtip = interp1(idtip_time, idtip, ip_prog_time)
            except mdsExceptions.MdsException:
                idtip = np.full(len(ip_prog_time), np.nan)

            try:
                istip, istip_time = params.mds_conn.get_data_with_dims(
                    r"\istip", tree_name="pcs_east"
                )
                istip = interp1(istip_time, istip, ip_prog_time)
            except mdsExceptions.MdsException:
                istip = np.full(len(ip_prog_time), np.nan)

            try:
                iutip, iutip_time = params.mds_conn.get_data_with_dims(
                    r"\iutip", tree_name="pcs_east"
                )
                iutip = interp1(iutip_time, iutip, ip_prog_time)
            except mdsExceptions.MdsException:
                iutip = np.full(len(ip_prog_time), np.nan)

            # Okay, now loop through all times, selecting the alternate isoflux target
            # Ip value for each time that isoflux control is enabled (i.e. for each
            # time that \sytps1 ~= 0).  (Only one of the 4 possible isoflux target Ip
            # signals is defined at each time.)
            tstart_isoflux_control = time_sytps1[sytps1_indices[0]]
            (indices,) = np.where(ip_prog_time >= tstart_isoflux_control)

            for it in range(len(indices)):
                if not np.isnan(ietip[it]):
                    ip_prog[it] = ietip[it]
                elif not np.isnan(idtip[it]):
                    ip_prog[it] = idtip[it]
                elif not np.isnan(istip[it]):
                    ip_prog[it] = istip[it]
                elif not np.isnan(iutip[it]):
                    ip_prog[it] = iutip[it]

            # End of block for dealing with isoflux control

        # For shots before year 2014, the LMTIPREF timebase needs to be shifted
        # by 17.0 ms
        if params.shot_id < 44432:
            ip_prog_time -= 0.0170

        # Calculate dipprog_dt
        dipprog_dt = np.gradient(ip_prog, ip_prog_time)

        # Interpolate all retrieved signals to the requested timebase
        ip = interp1(ip_time, ip, params.times)
        dip_dt = interp1(ip_time, dip_dt, params.time)
        ip_prog = interp1(ip_prog_time, ip_prog, params.times)
        dipprog_dt = interp1(ip_prog_time, dipprog_dt, params.times)

        # Calculate ip_error and ip_error_normalized
        ip_error = ip - ip_prog
        ip_error_normalized = ip_error / ip_prog

        output = {
            "ip": ip,
            "ip_prog": ip_prog,
            "ip_error": ip_error,
            "ip_error_normalized": ip_error_normalized,
            "dip_dt": dip_dt,
            "dipprog_dt": dipprog_dt,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["v_loop"],
        tokamak=Tokamak.EAST,
    )
    def get_v_loop(params: PhysicsMethodParams):
        """
        This routine gets the loop voltage from the EAST tree. The signal in the
        tree is derived by taking the time derivative of a flux loop near the
        inboard midplane.  Two possible signals are available, one digitised at a
        high rate (50 kHz), and the other sub-sampled down to 1 kHz.  This
        routine reads in the 1 kHz signal.  It linearly interpolates the loop
        voltage signal onto the specified timebase.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'v_loop' : array
                Calculated loop voltage [V].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_v_loop.m

        Original Author: Robert Granetz, Apr 2016

        Last major update: 11/19/24 by William Wei
        """
        v_loop = [np.nan]

        # Get "\vp1_s" signal from the EAST tree.  (This signal is a sub-sampled
        # version of "vp1".)
        v_loop, v_loop_time = params.mds_conn.get_data_with_dims(
            r"\vp1_s", tree_name="east"
        )

        # Interpolate the signal onto the requested timebase
        v_loop = interp1(v_loop_time, v_loop, params.times)

        return {"v_loop": v_loop}

    @staticmethod
    @physics_method(
        columns=[
            "zcur",
            "z_prog",
            "z_error",
            "zcur_lmsz",
            "z_error_lmsz",
            "zcur_lmsz_normalized",
            "z_error_lmsz_normalized",
            "aminor",
        ],
        tokamak=Tokamak.EAST,
    )
    def get_z_error(params: PhysicsMethodParams):
        """
        This script calculates Z_error = Z_cur - Z_programmed, or how much the
        actual vertical position differs from the requested position.  Two
        different methods are used to calculate Z_error, and both versions are
        returned by the routine.  The original method gets the z-centroid from
        the \zcur calculated by EFIT; the other uses the z-centroid from the
        \lmsz signal calculated by the plasma control system (PCS).  The
        normalized versions of the lmsz-derived signals are also returned.  They
        are normalized to the plasma minor radius.  All the signals are linearly
        interpolated onto the given timebase.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'zcur' : array
                Calculated z from EFIT [m].
            - 'z_prog' : array
                Programmed/requested/target z [m].
            - 'z_error' : array
                z_error = z_cur - z_prog [m].
            - 'zcur_lmsz' : array
                Calculated z from PCS [m].
            - 'z_error_lmsz' : array
                z_error_lmsz = zcur_lmsz - z_prog [m].
            - 'zcur_lmsz_normalized' : array
                zcur_lmsz / aminor
            - 'z_error_lmsz_normalized' : array
                z_error_lmsz_normalized = z_error_lmsz / aminor

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_Z_error.m

        Original Authors
        ----------------
        - Wang Bo, 2015/12/10
        - Alex Tinguely, 2015/09/09
        - Robert Granetz

        Last major update: 11/19/24 by William Wei
        """
        zcur = [np.nan]
        z_prog = [np.nan]
        z_error = [np.nan]
        zcur_lmsz = [np.nan]
        z_error_lmsz = [np.nan]
        zcur_lmsz_normalized = [np.nan]
        z_error_lmsz_normalized = [np.nan]
        aminor = [np.nan]

        # Read in the calculated zcur from EFIT
        zcur, zcur_time = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:zcur", tree_name="_efit_tree"
        )  # [A], [s]
        # Deal with rare bug
        zcur_time, unique_indices = np.unique(zcur_time, return_index=True)
        zcur = zcur[unique_indices]

        # Read in aminor from EFIT
        # TODO: use \aminor or \aout? -- MATLAB: \aminor
        aminor = params.mds_conn.get_data(r"\aminor", treename="_efit_tree")  # [m]
        aminor = aminor[unique_indices]

        # Next, get the programmed/requested/target Z from PCS, and the
        # calculated Z-centroid from PCS
        z_prog, z_prog_time = params.mds_conn.get_data_with_dims(
            r"\lmtzref/100", tree_name="pcs_east"
        )  # [m], [s] (Node says 'm' but it's wrong)
        zcur_lmsz, lmsz_time = params.mds_conn.get_data_with_dims(
            r"\lmsz", tree_name="pcs_east"
        )  # [m], [s]

        # Interpolate all retrieved signals to the requested timebase
        zcur = interp1(zcur_time, zcur, params.times)
        aminor = interp1(zcur_time, aminor, params.times)
        z_prog = interp1(z_prog_time, z_prog, params.times)
        zcur_lmsz = interp1(lmsz_time, zcur_lmsz, params.times)

        # Calculate both versions of z_error
        z_error = zcur - z_prog
        z_error_lmsz = zcur_lmsz - z_prog
        z_error_lmsz_normalized = z_error_lmsz / aminor

        output = {
            "zcur": zcur,
            "z_prog": z_prog,
            "z_error": z_error,
            "zcur_lmsz": zcur_lmsz,
            "z_error_lmsz": z_error_lmsz,
            "zcur_lmsz_normalized": zcur_lmsz_normalized,
            "z_error_lmsz_normalized": z_error_lmsz_normalized,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["ne", "greenwald_fraction", "dn_dt"],
        tokamak=Tokamak.EAST,
    )
    def get_density_parameters(params: PhysicsMethodParams):
        """
        This routine obtains the line-averaged density from the HCN vertical
        chord.  The node is called \DFSDEV in the PCS_EAST tree.  This signal is
        also used by the PCS for feedback control of the density.

        The plasma density is also measured by a POlarimeter/INTerferometer
        diagnostic (PO INT), which has 11 horizontal chords.  The midplane
        channel is \POINT_N6 in the EAST tree.  (This diagnostic uses an FIR
        laser at 432 micrometers.)  However, I found that the POINT density
        measurement is bad on a significant fraction of EAST shots, even in
        2017.

        The routine also calculates the Greenwald density using aminor coming
        from EFIT.  The time derivative of electron density is also obtained.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'ne' : array
                Density as measured by the polarimeter.interferometer [m^-3].
            - 'greenwald_fraction' : array
                greenwald_fraction = ne / nG [dimensionless].
            - 'dn/dt' : array
                dn_dt = d(ne)/dt [m^-3/s].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_density_parameters.m

        Original Author: Robert Granetz, Apr 2017

        Last major update: 11/19/24 by William Wei
        """
        ne = [np.nan]
        greenwald_fraction = [np.nan]
        dn_dt = [np.nan]

        # Get the density and calculate dn_dt
        ne, netime = params.mds_conn.get_data_with_dims(
            r"\dfsdev*1e19", tree_name="pcs_east"
        )  # [m^-3], [s]
        dn_dt = np.gradient(ne, netime)  # [m^-3/s]

        # Interpolate ne and dn_dt to the requested timebase
        ne = interp1(netime, ne, params.times)
        dn_dt = interp1(netime, dn_dt, params.times)

        # Calculate Greenwald density
        # TODO: use \aminor or \aout? -- MATLAB: \aout
        aminor, efittime = params.mds_conn.get_data_with_dims(
            r"\aout", tree_name="_efit_tree"
        )  # [m], [s]
        aminor = interp1(efittime, aminor, params.times)
        ip = EASTPhysicsMethods.get_ip_parameters(params)["ip"]  # [A]
        nG = 1e20 * (ip / 1e6) / (np.pi * aminor**2)  # [m^-3]
        greenwald_fraction = ne / nG

        output = {"ne": ne, "greenwald_fraction": greenwald_fraction, "dn_dt": dn_dt}
        return output

    @staticmethod
    @physics_method(
        columns=[
            "p_rad",
            "p_ecrh",
            "p_lh",
            "p_icrf",
            "p_nbi",
            "rad_input_frac",
            "rad_loss_frac",
            "p_input",
        ],
        tokamak=Tokamak.EAST,
    )
    def get_power(params: PhysicsMethodParams):
        """
        This function gets the input heating powers -- ohmic (p_OHM),electron
        cyclotron resonance heating (p_ECRH), neutral beam injection system (p_NBI)
        ion cyclotron (p_ICRF), and lower hybrid (p_LH) -- as well as the radiated
        output power (p_RAD).  If any of the auxiliary heating powers are not
        available (there was no ICRF or LH), then this function returns an array
        of zeros for them.  The ohmic heating power is obtained by calling a
        separate function, get_P_ohm.m.  If either the ohmic heating power or the
        radiated power is unavailable, then arrays of NaN (Not-a-Number) are
        returned for them, and for the radiated power fraction.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'p_rad' : array
                Radiated power [W].
            - 'p_ecrh' : array
                Electron cyclotron resonance heating power [W].
            - 'p_lh' : array
                Lower hybrid power [W].
            - 'p_icrf' : array
                Ion cyclotron resonance heating power [W].
            - 'p_nbi' : array
                Neutral beam injection power [W].
            - 'p_input' : array
                Total input power = p_ohm + p_lh + p_icrf + p_ecrh + p_nbi [W].
            - 'rad_input_frac' : array
                Radiated input fraction = p_rad / p_input [%].
            - 'rad_loss_frac' : array
                Radiated loss fraction = p_rad / (p_rad + dWmhd/dt) [%].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_power.m

        Original Author
        ----------------
        - Wang Bo, Dec 2015
        - Robert Granetz, Oct 2015
        - Alex Tinguely, Oct 2016

        Last major update: 11/20/24 by William Wei
        """
        p_rad = [np.nan]
        p_ecrh = [np.nan]
        p_lh = [np.nan]
        p_ohm = [np.nan]
        p_icrf = [np.nan]
        p_nbi = [np.nan]
        rad_input_frac = [np.nan]
        rad_loss_frac = [np.nan]
        p_input = [np.nan]

        def get_heating_power(addresses, tree):
            """
            Helper function to pack p_lh, p_icrf, and p_nbi computations
            """
            heating_power = np.empty(params.times.shape)
            for address in addresses:
                power_address, time_address = params.mds_conn.get_data_with_dims(
                    address, tree_name=tree
                )
                heating_power += interp1(
                    time_address,
                    power_address,
                    params.times,
                    kind="linear",
                    fill_value=0,
                )
            return heating_power

        # Get lower hybrid power

        # Note: the timebase for the LH power signal does not extend over the full
        # time span of the discharge.  Therefore, when interpolating the LH power
        # signal onto the "timebase" array, the LH signal has to be extrapolated
        # with zero values.  This is an option in the 'interp1' routine.  If the
        # extrapolation is not done, then the 'interp1' routine will assign NaN
        # (Not-a-Number) values for times outside the LH timebase, and the NaN's
        # will propagate into p_input, rad_input_frac, and rad_loss_frac, which is
        # not desirable.

        # TODO: check if the extrapolation method actually works
        lh_addresses = [
            r"\plhi1*1e3",  # LHW 2.45 GHz data
            r"\plhi2*1e3",  # LHW 4.66 GHz data
        ]
        p_lh = get_heating_power(lh_addresses, "east_1")  # [W]

        # Get ECRH power
        p_ecrh, ecrh_time = params.mds_conn.get_data_with_dims(
            r"\pecrh1i*1e3", tree_name="analysis"
        )  # [W], [s]
        (baseline_indices,) = np.where(ecrh_time < 0)
        if len(baseline_indices) > 0:
            p_ecrh_baseline = np.mean(p_ecrh[baseline_indices])
            p_ecrh -= p_ecrh_baseline
        p_ecrh = interp1(ecrh_time, p_ecrh, params.times, kind="linear", fill_value=0)

        # Get ICRF power
        # p_ICRF = (p_icrfii - p_icrfir) + (p_icrfbi - p_icrfbr)
        icrf_addresses = [
            r"\picrfii*1e3",
            r"\picrfir*-1e3",
            r"\picrfbi*1e3",
            r"\picrfbr*-1e3",
        ]
        p_icrf = get_heating_power(icrf_addresses, "icrf_east")  # [W]

        # Get NBI power
        nbi_addresses = [
            r"\pnbi1lsource*1e3",
            r"\pnbi1rsource*1e3",
            r"\pnbi2lsource*1e3",
            r"\pnbi2rsource*1e3",
        ]
        p_nbi = get_heating_power(nbi_addresses, "nbi_east")  # [W]

        # Get ohmic power
        p_ohm = EASTPhysicsMethods.get_p_ohm(params)["p_ohm"]  # [W]

        # Get radiated power
        # tree = 'prad_east'
        # TODO: find and implement `prad_bulk_xuv2014_2016`
        # p_rad, rad_time = prad_bulk_xuv2014_2016(params.shot_id)    # Original from Duan Yanming, and modified by RSG
        # p_rad = interp1(rad_time, p_rad, params.times, kind="linear", fill_value=0)
        p_rad = [np.nan] * len(params.times)

        # Get Wmhd and calculate dWmhd_dt
        wmhd, efittime = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:wplasm", tree_name="_efit_tree"
        )  # [W], [s]
        dwmhd_dt = np.gradient(wmhd, efittime)
        dwmhd_dt = interp1(efittime, dwmhd_dt, params.times)

        # Calculate p_input, rad_input_frac, and rad_loss_frac
        p_input = p_ohm + p_lh + p_icrf + p_ecrh + p_nbi
        rad_input_frac = (
            p_rad / p_input
        )  # TODO: div/0 error? Use with np.errorstate(...)
        rad_loss_frac = p_rad / (p_input - dwmhd_dt)

        output = {
            "p_rad": p_rad,
            "p_ecrh": p_ecrh,
            "p_lh": p_lh,
            "p_ohm": p_ohm,
            "p_icrf": p_icrf,
            "p_nbi": p_nbi,
            "rad_input_frac": rad_input_frac,
            "rad_loss_frac": rad_loss_frac,
            "p_input": p_input,
        }
        return output

    @staticmethod
    @physics_method(
        columns=["p_ohm"],
        tokamak=Tokamak.EAST,
    )
    def get_p_ohm(params: PhysicsMethodParams):
        """
        This script calculates the ohmic power, p_ohm. We use the following
        expression to calculate P_ohm:

        P_ohm = Ip * V_resistive

        where: V_resistive = V_loop - V_inductive = V_loop - L * dIp/dt
        and L = L_internal = mu0 * R0 * li/2

        - pcs_east node '\pcvloop'        is used for V_loop,
        - efit18   node '\efit_aeqdsk:li' is used for li,
        - pcs_east node '\pcrl01'         is used for Ip

        If the EFIT data or the magnetics Ip data is not available, then P_ohm is
        returned as NaN.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'p_ohm' : array
                Ohmic power [W].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_P_ohm.m

        Original Authors
        ----------------
        Wang Bo, Dec 2015
        Alex Tinguely, Sep 2015
        Robert Granetz, Oct 2015 -- Jun 2016

        Last major update: 2014/11/21 by William Wei
        """
        p_ohm = [np.nan]

        # Get raw signals
        vloop, vloop_time = params.mds_conn.get_data_with_dims(
            r"\pcvloop", tree_name="pcs_east"
        )  # [V]
        li, li_time = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:li", tree_name="_efit_tree"
        )  # [H]
        # Fetch raw ip signal to calculate dip_dt and apply smoothing
        ip, ip_time = params.mds_conn.get_data_with_dims(
            r"\pcrl01", tree_name="pcs_east"
        )  # [A]
        sign_ip = np.sign(sum(ip))
        dipdt = np.gradient(ip, ip_time)
        dipdt_smoothed = smooth(dipdt, 11)  # Use 11-point boxcar smoothing

        # Interpolate fetched signals to the requested timebase
        for signal, signal_time in zip(
            [vloop, li, ip, dipdt_smoothed], [vloop_time, li_time, ip_time, ip_time]
        ):
            signal = interp1(
                signal_time, signal, params.times, kind="linear", bounds_error=0
            )

        # Calculate p_ohm
        # TODO: check DIV/0 error
        # TODO: replace 1.85 with fetched r0 signal
        inductance = 4e-7 * np.pi * 1.85 * li / 2  # For EAST use R0 = 1.85 m
        v_inductive = -inductance * dipdt_smoothed
        v_resistive = vloop - v_inductive
        (v_resistive_indices,) = np.where(v_resistive * sign_ip < 0)
        v_resistive[v_resistive_indices] = 0
        p_phm = ip * v_resistive

        return {"p_ohm": p_ohm}

    @staticmethod
    @physics_method(
        columns=[
            "n_equal_1_mode",
            "n_equal_1_phase",
            "e_equal_1_normalized",
            "rmp_n_equal_1",
            "rmp_n_equal_1_phase",
            "btor",
        ],
        tokamak=Tokamak.EAST,
    )
    def get_n_equal_1_data(params: PhysicsMethodParams):
        """
        This function computes the amplitude and phase of the n=1 Fourier
        component of the net saddle signals (total saddle signals minus the
        calculated pickup from the RMP coils) and interpolates the signals onto
        the specified timebase.  It also computes the amplitude and phase of the
        n=1 Fourier component of the calculated pickup from the RMP coils into
        the database.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'n_equal_1_mode' : array
                Amplitude of n=1 Fourier component of saddle signals after
                subtracting the calculated pickup from the RMP coil currents [T].
            - 'n_equal_1_phase' : array
                Toroidal phase angle of above [rad].
            - 'n_equal_1_normalized' : array
                'n_equal_1_mode' normalized to btor.
            'rmp_n_equal_1' : array
                Amplitude of the n=1 Fourier component of the calculated pickup
                of the RMP coils on the saddle signals [T].
            'rmp_n_equal_1_phase' : array
                toroidal phase angle of above [rad].
            'btor' : array
                Toroidal magnetic field [T].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_n_equal_1_data.m

        Original Authors
        ----------------
        Robert Granetz, Apr 2017

        Last major update: 2014/11/21 by William Wei
        """
        n_equal_1_mode = [np.nan]
        n_equal_1_phase = [np.nan]
        n_equal_1_normalized = [np.nan]
        rmp_n_equal_1 = [np.nan]
        rmp_n_equal_1_phase = [np.nan]
        btor = [np.nan]

        # Read in the saddle sensor data and rmp coil currents
        # TODO: implement get_rmp_and_saddle_signals
        # rmptime, rmp, saddletime, saddle = get_rmp_and_saddle_signals(params.shot_id)

        output = {
            "n_equal_1_mode": n_equal_1_mode,
            "n_equal_1_phase": n_equal_1_phase,
            "n_equal_1_normalized": n_equal_1_normalized,
            "rmp_n_equal_1": rmp_n_equal_1,
            "rmp_n_equal_1_phase": rmp_n_equal_1_phase,
            "btor": btor,
        }

        return output

    @staticmethod
    @physics_method(
        columns=["upper_gap", "lower_gap"],
        tokamak=Tokamak.EAST,
    )
    def get_efit_gaps(params: PhysicsMethodParams):
        """
        This script calculates upper and lower gaps from EFIT information
        on plasma boundary and first wall geometry. This is needed because
        the EAST version of EFIT does not output the upper and lower gaps
        directly.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'upper_gap' : array
                Upper gap [m].
            - 'lower_gap' : array
                Lower gap [m].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_EFIT_gaps.m

        Original Authors
        ----------------
        Robert Granetz
        Jiaxiang Zhu

        Last major update: 2014/11/22 by William Wei
        """
        upper_gap = [np.nan]
        lower_gap = [np.nan]

        # TODO: verify all of the reshape and tile functions!
        # Get plasma boundary data
        data, efittime = params.mds_conn.get_data_with_dims(
            r"\top.results.geqdsk:bdry", tree_name="_efit_tree"
        )
        # TODO: Check the actual shape in MATLAB script
        xcoords = np.reshape(data[0, :, :], (-1, len(efittime)))
        ycoords = np.reshape(data[1, :, :], (-1, len(efittime)))

        # Get first wall geometry data
        xfirstwall = params.mds_conn.get_data(
            r"\top.results.geqdsk:xlim", tree_name="_efit_tree"
        )
        yfirstwall = params.mds_conn.get_data(
            r"\top.results.geqdsk:ylim", tree_name="_efit_tree"
        )
        seed = np.ones((len(xcoords), 1))
        xfirstwall = np.reshape(xfirstwall, (-1, 1))
        yfirstwall = np.reshape(yfirstwall, (-1, 1))
        xfirstwall_mat = np.tile(xfirstwall, (1, len(efittime)))
        yfirstwall_mat = np.tile(yfirstwall, (1, len(efittime)))

        # Calculate upper & lower gaps
        (index_upperwall,) = np.where(yfirstwall > 0.6)
        (index_lowerwall,) = np.where(yfirstwall < -0.6)

        xupperwall = xfirstwall_mat[index_upperwall, :]
        xupperwall_mat = np.reshape(
            np.kron(xupperwall, seed),
            (
                len(
                    seed,
                ),
                -1,
                len(efittime),
            ),
        )
        xlowerwall = xfirstwall_mat[index_lowerwall, :]
        xlowerwall_mat = np.reshape(
            np.kron(xlowerwall, seed),
            (
                len(
                    seed,
                ),
                -1,
                len(efittime),
            ),
        )

        yupperwall = yfirstwall_mat[index_upperwall, :]
        yupperwall_mat = np.reshape(
            np.kron(yupperwall, seed),
            (
                len(
                    seed,
                ),
                -1,
                len(efittime),
            ),
        )
        ylowerwall = yfirstwall_mat[index_lowerwall, :]
        ylowerwall_mat = np.reshape(
            np.kron(ylowerwall, seed),
            (
                len(
                    seed,
                ),
                -1,
                len(efittime),
            ),
        )

        xupperplasma_mat = np.reshape(
            np.tile(xcoords, (len(xupperwall), 1))(-1, len(xupperwall), len(efittime))
        )
        yupperplasma_mat = np.reshape(
            np.tile(ycoords, (len(yupperwall), 1))(-1, len(yupperwall), len(efittime))
        )
        xlowerplasma_mat = np.reshape(
            np.tile(xcoords, (len(xlowerwall), 1))(-1, len(xlowerwall), len(efittime))
        )
        ylowerplasma_mat = np.reshape(
            np.tile(ycoords, (len(ylowerwall), 1))(-1, len(ylowerwall), len(efittime))
        )

        uppergap_mat = (
            ((xupperplasma_mat - xupperwall_mat) ** 2)
            + ((yupperplasma_mat - yupperwall_mat) ** 2)
        ) ** 0.5
        lowergap_mat = (
            ((xlowerplasma_mat - xlowerwall_mat) ** 2)
            + ((ylowerplasma_mat - ylowerwall_mat) ** 2)
        ) ** 0.5
        upper_gap = np.reshape(min(min(uppergap_mat)), (-1, 1))
        lower_gap = np.reshape(min(min(lowergap_mat)), (-1, 1))

        # Interpolate to the requested timebase
        # TODO: do this after clean up all the above codes!

        return {"upper_gap": upper_gap, "lower_gap": lower_gap}

    @staticmethod
    @physics_method(
        columns=["pupper_gap", "plower_gap"],
        tokamak=Tokamak.EAST,
    )
    def get_pefit_gaps(params: PhysicsMethodParams):
        """
        This script calculates upper and lower gaps from P-EFIT information
        on plasma boundary and first wall geometry. This is needed because
        the EAST version of EFIT does not output the upper and lower gaps
        directly.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'pupper_gap' : array
                Upper gap computed from P-EFIT [m].
            - 'plower_gap' : array
                Lower gap computed from P-EFIT [m].

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_pEFIT_gaps.m

        Original Authors
        ----------------
        Robert Granetz
        Jiaxiang Zhu

        Last major update: 2014/11/22 by William Wei
        """
        upper_gap = [np.nan]
        lower_gap = [np.nan]

        return {"pupper_gap": upper_gap, "plower_gap": lower_gap}

    @staticmethod
    @physics_method(
        columns=["ip_error_rt", "q95_rt", "beta_p_rt", "li_rt", "wmhd_rt"],
        tokamak=Tokamak.EAST,
    )
    def get_pcs_parameters(params: PhysicsMethodParams):
        """
        This function gets many of the real time signals that are actually used
        in the plasma control system (PCS) on the EAST tokamak.

        For development of a disruption prediction algorithm that we want to run
        in real time in the PCS, it is better to train on a database of these real
        time signals than on the processed signals in the analysis trees and east
        trees.

        Note: since 2018 the EFIT-derived signals used in the PCS are calculated
        by P-EFIT, which is different than RT-EFIT.  This information came from
        QP Yuan <qpyuan@ipp.ac.cn>

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'ip_error_rt' : array
                error between actual and pre-programmed plasma currents
            - 'q95_rt' : array
                q95 calculated by RT-EFIT for PCS (P-EFIT since 2018)
            - 'beta_p_rt' : array
                beta_p calculated by RT-EFIT for PCS (P-EFIT since 2018)
            - 'li_rt' : array
                li calculated by RT-EFIT for PCS (P-EFIT since 2018)
            - 'wmhd_rt' : array
                Wmhd calculated by RT-EFIT for PCS (P-EFIT since 2018)

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_pcs.m

        Original Author
        ----------------
        Robert Granetz, Dec 2018

        Last major update: 2014/11/22 by William Wei
        """
        output = dict()

        # Get ip_error_rt, beta_p_rt, li_rt, and wmhd_rt
        signals = {
            "ip_error_rt": r"\lmeip",
            "beta_p_rt": r"\pfsbetap",
            "li_rt": r"\pfsli",
            "wmhd_rt": r"\pfswmhd",
        }
        for name, address in signals.items():
            signal, timearray = params.mds_conn.get_data_with_dims(
                address, tree_name="pcs_east"
            )
            signal = interp1(
                timearray, signal, params.times, kind="linear", bounds_error=0
            )
            output[name] = signal

        # Get q95_rt
        q95_rt, q95_rt_time = params.mds_conn.get_data_with_dims(
            r"\q95", tree_name="pefitrt_east"
        )
        # Deal with bug
        q95_rt_time, unique_indices = np.unique(q95_rt_time, return_index=True)
        q95_rt = q95_rt[unique_indices]
        q95_rt = interp1(
            q95_rt_time, q95_rt, params.times, kind="linear", bounds_error=0
        )
        output["q95_rt"] = q95_rt

        return output

    @staticmethod
    @physics_method(
        columns=[
            "p_rad_rt",
            "p_lh_rt",
            "p_nbi_rt",
        ],
        tokamak=Tokamak.EAST,
    )
    def get_pcs_power(params: PhysicsMethodParams):
        """
        This function gets the real time power signals that are actually used
        in the plasma control system (PCS) on the EAST tokamak.

        For development of a disruption prediction algorithm that we want to run
        in real time in the PCS, it is better to train on a database of these real
        time signals than on the processed signals in the analysis trees and east
        trees.

        Note: since 2018 the EFIT-derived signals used in the PCS are calculated
        by P-EFIT, which is different than RT-EFIT.  This information came from
        QP Yuan <qpyuan@ipp.ac.cn>

        As of Nov 2024, only p_rad_rt, p_lh_rt, and p_nbi_rt have been implemented.
        All other PCS signals weren't implemented in original MATLAB script.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'p_rad_rt' : array
                radiated power [W]
            - [UNAVAILABLE] 'p_ecrh_rt' : array
                electron cyclotron resonance heating power [W]
            - 'p_lh_RT' : array
                lower hybrid power [W]
            - [UNAVAILABLE] 'p_oh_rt' : array
                ohmic power [W]
            - [UNAVAILABLE] 'p_icrf_rt' : array
                ion cyclotron power [W]
            - 'p_nbi_RT' : array
                neutral beam injection power [W]
            - [UNAVAILABLE] 'rad_input_frac' : array
                rad_input_frac = p_rad/p_input [%]
                               = p_rad/(p_ohm + p_lh + p_icrh + p_ecrh + p_nbi)
                               = ratio of radiated power to total input power
            - [UNAVAILABLE] 'rad_loss_frac' : array
                rad_loss_frac = p_rad/p_loss
                              = p_rad/(p_rad + p_cond + p_conv)
                              = p_rad/(p_input - dWmhd/dt)
                              = ratio of radiated power to total loss power

        References
        -------
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_pcs.m

        Original Author
        ----------------
        Robert Granetz, Dec 2018

        Last major update: 2014/11/22 by William Wei
        """
        # Get p_rad_rt
        p_rad_rt, timearray = params.mds_conn.get_data_with_dims(
            r"\pcprad", tree_name="pcs_east"
        )
        p_rad_rt = interp1(
            timearray, p_rad_rt, params.times, kind="linear", bounds_error=0
        )

        # TODO: Verify the outputs, then modify get_heating_power() to make it compatible with this
        # Get p_nbi_rt
        nbi_addresses = {
            "i_nbil_rt": r"\pcnbi1li",
            "v_nbil_rt": r"\pcnbi1lv",
            "i_nbir_rt": r"\pcnbi1ri",
            "v_nbir_rt": r"\pcnbi1rv",
        }
        nbi_signals = dict()
        for name, address in nbi_addresses.items():
            signal, timearray = params.mds_conn.get_data_with_dims(
                address, tree_name="pefitrt_east"
            )
            signal = interp1(
                timearray, signal, params.times, kind="linear", bounds_error=0
            )
            nbi_signals[name] = signal
        p_nbi_rt = (
            nbi_signals["i_nbil_rt"] * nbi_signals["v_nbil_rt"]
            + nbi_signals["i_nbir_rt"] * nbi_signals["v_nbir_rt"]
        )

        # Get p_lh_rt
        lh_addresses = {
            "p_lh_46_inj_rt": r"\pcplhi",
            "p_lh_46_ref_rt": r"\pcplhr",
            "p_lh_245_inj_rt": r"\pcplhi2",
            "p_lh_245_ref_rt": r"\pcplhr2",
        }
        lh_signals = dict()
        for name, address in lh_addresses.items():
            signal, timearray = params.mds_conn.get_data_with_dims(
                address, tree_name="pefitrt_east"
            )
            signal = interp1(
                timearray, signal, params.times, kind="linear", bounds_error=0
            )
            lh_signals[name] = signal
        p_lh_rt = (
            lh_signals["p_lh_46_inj_rt"]
            - lh_signals["p_lh_46_ref_rt"]
            + lh_signals["p_lh_245_inj_rt"]
            - lh_signals["p_lh_245_ref_rt"]
        )

        # Q.P. Yuan:
        # There is no signals from ICRF or ECRH connected to PCS. And I checked the
        # signals for NBI, there is no effective value, just noise. We will make
        # the NBI signals available in this EAST campaign. The PLHI2 and PLHR2 for
        # 2.45G LHW system have signals but are not calibrated yet.

        # Get p_ohm
        # "Note, I'm not bothering to calculate a real time version" -- Whoever who wrote this
        # p_ohm = EASTPhysicsMethods.get_p_ohm(params)['p_ohm']

        # Compute total power and radiation fractions
        # p_input_rt = p_ohm + p_lh_rt + p_icrf_rt + p_ecrh_rt + p_nbi_rt
        # rad_input_frac_rt = p_rad_rt / p_input_rt
        # rad_loss_frac_rt = p_rad_rt / (p_input_rt - dwmhd_dt)

        output = {
            "p_rad_rt": p_rad_rt,
            "p_lh_rt": p_lh_rt,
            "p_nbi_rt": p_nbi_rt,
        }
        return output
