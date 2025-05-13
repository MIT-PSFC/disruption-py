#!/usr/bin/env python3

"""
Module for retrieving and calculating data for DIII-D physics methods.
"""

import numpy as np
import scipy
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1, matlab_smooth
from disruption_py.machine.east.efit import EastEfitMethods
from disruption_py.machine.east.util import EastUtilMethods
from disruption_py.machine.tokamak import Tokamak


class EastPhysicsMethods:
    """
    A class to retrieve and calculate physics-related data for EAST.
    """

    @staticmethod
    @physics_method(columns=["time_until_disrupt"], tokamak=Tokamak.EAST)
    def get_time_until_disrupt(params: PhysicsMethodParams):
        """
        Calculate the time until disruption.

        Currently, the disruption time is queried from the `DISRUPTIONS` table
        in the SQL database of each machine. These disruption times were calculated
        using Robert Granetz's routine.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the disruption information and times.

        Returns
        -------
        dict
            A dictionary with a single key `time_until_disrupt`.

        References
        -------
        - issues: #[223](https://github.com/MIT-PSFC/disruption-py/issues/223)
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
        This routine retrieves the plasma current signals. It then calculates the error
        between the programmed and the actual plasma current, and the time derivatives
        of the actual and the programmed plasma current.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `ip`: Measured plasma current.
            - `ip_prog`: Programmed plasma current.
            - `ip_error`: Error between the actual and programmed plasma currents.
            - `ip_error_normalized`: `ip_error` normalized to `ip_prog`.
            - `dip_dt`: Time derivative of the measured plasma current.
            - `dipprog_dt`: Time derivative of the programmed plasma current.

        References
        -------
        - original source: [get_Ip_parameters.m](https://github.com/MIT-PSFC/disruption-py
        /blob/matlab/EAST/get_Ip_parameters.m)
        """
        ip = [np.nan]
        ip_prog = [np.nan]
        ip_error = [np.nan]
        ip_error_normalized = [np.nan]
        dip_dt = [np.nan]
        dipprog_dt = [np.nan]

        ip, ip_time = EastUtilMethods.retrieve_ip(params.mds_conn, params.shot_id)
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
        dip_dt = interp1(ip_time, dip_dt, params.times)
        ip_prog = interp1(ip_prog_time, ip_prog, params.times)
        dipprog_dt = interp1(ip_prog_time, dipprog_dt, params.times)

        # Calculate ip_error and ip_error_normalized
        ip_error = ip - ip_prog
        ip_error_normalized = ip_error / ip_prog

        return {
            "ip": ip,
            "ip_prog": ip_prog,
            "ip_error": ip_error,
            "ip_error_normalized": ip_error_normalized,
            "dip_dt": dip_dt,
            "dipprog_dt": dipprog_dt,
        }

    @staticmethod
    @physics_method(columns=["v_loop"], tokamak=Tokamak.EAST)
    def get_v_loop(params: PhysicsMethodParams):
        r"""
        Get the loop voltage signal.

        By default, this routine gets the loop voltage from \vp1_s from the east tree.
        The signal in the tree is derived by taking the time derivative of a flux loop
        near the inboard midplane.  Two possible signals are available, one digitised at a
        high rate (50 kHz), and the other sub-sampled down to 1 kHz.  This
        routine reads in the 1 kHz signal.

        If \vp1_s isn't available, the method will fall back to using \pcvloop from the
        pcs_east tree instead.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the loop voltage (`v_loop`).

        References
        -------
        - original source: [get_v_loop.m](https://github.com/MIT-PSFC/disruption-py
        /blob/matlab/EAST/get_v_loop.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411),
        #[451](https://github.com/MIT-PSFC/disruption-py/pull/451)
        """
        v_loop = [np.nan]

        # Get "\vp1_s" signal from the EAST tree.  (This signal is a sub-sampled
        # version of "vp1".)
        try:
            v_loop, v_loop_time = params.mds_conn.get_data_with_dims(
                r"\vp1_s", tree_name="east"
            )
        except mdsExceptions.MdsException:
            params.logger.verbose(
                r"v_loop: Failed to get \vp1_s data. Use \pcvloop from pcs_east instead."
            )
            v_loop, v_loop_time = params.mds_conn.get_data_with_dims(
                r"\pcvloop", tree_name="pcs_east"
            )  # [V]

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
        ],
        tokamak=Tokamak.EAST,
    )
    def get_z_error(params: PhysicsMethodParams):
        r"""
        Get the programmed and measured vertical position of the plasma current
        centroid, then calculate the control error.

        This routine uses two methods to compute z_error. The original method
        gets the z-centroid from the \zcur calculated by EFIT; the other uses
        the z-centroid from the \lmsz signal calculated by the plasma control
        system (PCS).

        This routine also returns the normalized versions of lmsz-derived signals.
        These signals are normalized to the plasma minor radius.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `zcur`: Calculated z from EFIT [m].
            - `z_prog`: Programmed/requested/target z [m].
            - `z_error`: z_error = z_cur - z_prog [m].
            - `zcur_lmsz`: Calculated z from PCS [m].
            - `z_error_lmsz`: zcur_lmsz - z_prog [m].
            - `zcur_lmsz_normalized`: zcur_lmsz / aminor.
            - `z_error_lmsz_normalized`: z_error_lmsz / aminor.

        References
        -------
        - original source: [get_Z_error.m](https://github.com/MIT-PSFC/disruption-py/
        blob/matlab/EAST/get_Z_error.m)
        """
        # Read in the calculated zcur from EFIT
        zcur, zcur_time = params.mds_conn.get_data_with_dims(
            r"\efit_aeqdsk:zcur", tree_name="_efit_tree"
        )  # [A], [s]
        # Deal with rare bug
        zcur_time, unique_indices = np.unique(zcur_time, return_index=True)
        zcur = zcur[unique_indices]

        # Read in aminor from EFIT
        # TODO: use \aminor or \aout? -- MATLAB: \aminor
        aminor = params.mds_conn.get_data(r"\aminor", tree_name="_efit_tree")  # [m]
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

        # Calculate zcur_lmsz_normalized
        zcur_lmsz_normalized = zcur_lmsz / aminor

        # Calculate both versions of z_error and z_error_lmsz_normalized
        z_error = zcur - z_prog
        z_error_lmsz = zcur_lmsz - z_prog
        z_error_lmsz_normalized = z_error_lmsz / aminor

        return {
            "zcur": zcur,
            "z_prog": z_prog,
            "z_error": z_error,
            "zcur_lmsz": zcur_lmsz,
            "z_error_lmsz": z_error_lmsz,
            "zcur_lmsz_normalized": zcur_lmsz_normalized,
            "z_error_lmsz_normalized": z_error_lmsz_normalized,
        }

    @staticmethod
    @physics_method(
        columns=["n_e", "greenwald_fraction", "dn_dt"], tokamak=Tokamak.EAST
    )
    def get_density_parameters(params: PhysicsMethodParams):
        r"""
        Calculate electron density, its time derivative, and the Greenwald fraction.

        The Greenwald fraction is the ratio of the measured electron density $n_e$ and
        the Greenwald density limit $n_G$ defined as [^1]:

        $$
        n_G = \frac{I_p}{\pi a^2}
        $$

        where $n_G$ is given in $10^{20} m^{-3}$ and $I_p$ is in MA.

        The line-averaged electron density is obtained from the HCN vertical chord.
        This signal is also used by the PCS for feedback control of the density.

        [^1]: https://wiki.fusion.ciemat.es/wiki/Greenwald_limit

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing electron density (`n_e`), its time derivative (`dn_dt`),
            and the Greenwald fraction (`greenwald_fraction`).

        References
        -------
        - original source: [get_density_parameters.m](https://github.com/MIT-PSFC/disruption-py
        /blob/matlab/EAST/get_density_parameters.m)
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
            r"\efit_aeqdsk:aout", tree_name="_efit_tree"
        )  # [m], [s]
        aminor = interp1(efittime, aminor, params.times)
        ip = EastPhysicsMethods.get_ip_parameters(params)["ip"]  # [A]
        n_g = 1e20 * (ip / 1e6) / (np.pi * aminor**2)  # [m^-3]
        greenwald_fraction = ne / n_g

        return {"n_e": ne, "greenwald_fraction": greenwald_fraction, "dn_dt": dn_dt}

    @staticmethod
    def _get_raw_axuv_data(params: PhysicsMethodParams):
        """
        Get the raw (uncalibrated) data from the AXUV arrays for calculating
        the total radiated power and prad_peaking. The AXUV diagnostic contains
        4 arrays of 16 channels each.

        Note that the original implementations in Prad_bulk_xuv2014_2016.m and
        get_prad_peaking.m are slightly different in where the calibration factors
        are applied. This difference isn't explained in the scripts so I keep
        these functions different.
        """
        # Get XUV data
        (xuvtime,) = params.mds_conn.get_dims(r"\pxuv1", tree_name="east_1")  # [s]
        # There are 64 AXUV chords, arranged in 4 arrays of 16 channels each
        xuv = np.full((len(xuvtime), 64), np.nan)

        dt = xuvtime[1] - xuvtime[0]
        smoothing_time = 1e-3
        smoothing_window = int(max([round(smoothing_time / dt, 1)]))

        for iarray in range(4):
            for ichan in range(16):
                ichord = 16 * iarray + ichan
                signal = params.mds_conn.get_data(
                    r"\pxuv" + str(ichord + 1), tree_name="east_1"
                )
                # Subtract baseline
                signal = signal - np.mean(signal[:100])
                # TODO: change this to causal smoothing
                xuv[:, ichord] = (
                    matlab_smooth(signal, smoothing_window) * 1e3
                )  # [kW] -> [W]
        return xuv, xuvtime

    @staticmethod
    @physics_method(columns=["p_rad"], tokamak=Tokamak.EAST)
    def get_radiated_power(params: PhysicsMethodParams):
        """
        Calculate total radiated power in bulk plasma from the AXUV arrays.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the total radiated power (`p_rad`).

        References
        -------
        - original source: Prad_bulk_xuv2014_2016.m (Currently not available in the repository)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411)
        """
        # Get the raw AXUV data
        try:
            xuv, xuvtime = EastPhysicsMethods._get_raw_axuv_data(params)  # [W], [s]
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get raw AXUS data")
            return {"p_rad": [np.nan]}

        # Get calibration factors
        calib_factors = EastUtilMethods.get_axuv_calib_factors()
        # Apply calibration factors Fac1 to Fac4
        for i in range(64):
            xuv[:, i] *= (
                calib_factors["Fac1"][i % 16]
                * calib_factors["Fac2"][i // 16][i % 16]
                * calib_factors["Fac3"][i // 16]
                * calib_factors["Fac4"]
            )
        # Correction for bad channels (from Duan Yanmin's program)
        # TODO: why not [35] = ([34]+[36])/2 when [35] is the bad channel?
        xuv[:, 11] = 0.5 * (xuv[:, 10] + xuv[:, 12])
        xuv[:, 35] = 0.5 * (xuv[:, 34] + xuv[:, 35])
        # Apply geometric calibrations and Fac5
        for i in range(64):
            xuv[:, i] *= (
                2
                * np.pi
                * calib_factors["Maj_R"]
                * calib_factors["Del_r"][i]
                * calib_factors["Fac5"]
            )

        # Subtract divertor measurements and overlapping measurments
        core_indices = (
            list(range(7, 15))
            + list(range(20, 32))
            + list(range(32, 48))
            + list(range(53, 59))
        )
        p_rad = np.nansum(xuv[:, core_indices], axis=1)
        # Interpolate to requested time base
        p_rad = interp1(xuvtime, p_rad, params.times, kind="linear", fill_value=0)
        return {"p_rad": p_rad}

    @staticmethod
    @physics_method(
        columns=[
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
        This routine gets the auxiliary heating powers: electron cyclotron
        resonance heating (p_ecrh), neutral beam injection system (p_nbi)
        ion cyclotron (p_icrf)), and lower hybrid (p_lh). If any of the auxiliary
        heating powers were not available during a shot, then it returns an array
        of zeros for them.

        This routine also calculates the total input power, the radiated input
        fraction, and the radiated loss function by calling get_p_ohm and
        get_radiated_power.

        Note that currently we assume that the MDSplus trees always have the corresponding
        tree nodes for each of the powers, even if a heating source/diagnostic was turned
        off in that shot. If the routine fails to find a tree or a tree node,
        then it assumes that tree is broken and returns the corresponding power as an array
        of nans. The routine will also skip calculating p_input or the radiated fractions if
        any of the powers used in the computation is an array of nans.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `p_ecrh`: Electron cyclotron resonance heating power [W].
            - `p_lh`: Lower hybrid power [W].
            - `p_icrf`: Ion cyclotron resonance heating power [W].
            - `p_nbi`: Neutral beam injection power [W].
            - `p_input`: Total input power = p_ohm + p_lh + p_icrf + p_ecrh + p_nbi [W].
            - `rad_input_frac: Radiated input fraction = p_rad / p_input [%].
            - `rad_loss_frac`: Radiated loss fraction = p_rad / (p_rad + dWmhd/dt) [%].

        References
        -------
        - original source: [get_power.m](https://github.com/MIT-PSFC/disruption-py/blob/
        matlab/EAST/get_power.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411), #[451](https:
        //github.com/MIT-PSFC/disruption-py/pull/451)
        """
        p_ecrh = [np.nan]
        p_lh = [np.nan]
        p_icrf = [np.nan]
        p_nbi = [np.nan]
        rad_input_frac = [np.nan]
        rad_loss_frac = [np.nan]
        p_input = [np.nan]

        def get_heating_power(nodes, tree):
            """
            Helper function to pack p_lh, p_icrf, and p_nbi computations
            """
            heating_power = np.zeros(params.times.shape)
            for node in nodes:
                power_node, time_node = params.mds_conn.get_data_with_dims(
                    node, tree_name=tree
                )
                heating_power += interp1(
                    time_node,
                    power_node,
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
        lh_nodes = [
            r"\plhi1*1e3",  # LHW 2.45 GHz data
            r"\plhi2*1e3",  # LHW 4.66 GHz data
        ]
        try:
            p_lh = get_heating_power(lh_nodes, "east_1")  # [W]
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get LH heating power")
            p_lh = [np.nan]

        # Get ECRH power
        try:
            p_ecrh, ecrh_time = params.mds_conn.get_data_with_dims(
                r"\pecrh1i*1e3", tree_name="analysis"
            )  # [W], [s]
            (baseline_indices,) = np.where(ecrh_time < 0)
            if len(baseline_indices) > 0:
                p_ecrh_baseline = np.mean(p_ecrh[baseline_indices])
                p_ecrh -= p_ecrh_baseline
            p_ecrh = interp1(
                ecrh_time, p_ecrh, params.times, kind="linear", fill_value=0
            )
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get ECRH heating power")
            p_ecrh = [np.nan]

        # Get ICRF power
        # p_ICRF = (p_icrfii - p_icrfir) + (p_icrfbi - p_icrfbr)
        icrf_nodes = [
            r"\picrfii*1e3",
            r"\picrfir*-1e3",
            r"\picrfbi*1e3",
            r"\picrfbr*-1e3",
        ]
        try:
            p_icrf = get_heating_power(icrf_nodes, "icrf_east")  # [W]
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get ICRF heating power")
            p_icrf = [np.nan]

        # Get NBI power
        nbi_nodes = [
            r"\pnbi1lsource*1e3",
            r"\pnbi1rsource*1e3",
            r"\pnbi2lsource*1e3",
            r"\pnbi2rsource*1e3",
        ]
        try:
            p_nbi = get_heating_power(nbi_nodes, "nbi_east")  # [W]
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get NBI heating power")
            p_nbi = [np.nan]

        # Get ohmic power
        p_ohm = EastPhysicsMethods.get_p_ohm(params)["p_oh"]  # [W]

        # Get radiated power
        p_rad = EastPhysicsMethods.get_radiated_power(params)["p_rad"]  # [W]

        # Calculate p_input
        if ~(
            np.isnan(p_ohm).all()
            or np.isnan(p_lh).all()
            or np.isnan(p_icrf).all()
            or np.isnan(p_ecrh).all()
            or np.isnan(p_nbi).all()
        ):
            p_input = p_ohm + p_lh + p_icrf + p_ecrh + p_nbi
        else:
            p_input = [np.nan]

        # Get Wmhd, calculate dWmhd_dt, and calculate p_loss
        try:
            wmhd, efittime = params.mds_conn.get_data_with_dims(
                r"\efit_aeqdsk:wmhd", tree_name="_efit_tree"
            )  # [W], [s]
            dwmhd_dt = np.gradient(wmhd, efittime)
            dwmhd_dt = interp1(efittime, dwmhd_dt, params.times)
            if not np.isnan(p_input).all():
                p_loss = p_input - dwmhd_dt
            else:
                p_loss = [np.nan]
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get proper signals to compute p_loss")
            p_loss = [np.nan]

        rad_input_frac = np.full(len(params.times), np.nan)
        rad_loss_frac = np.full(len(params.times), np.nan)
        if ~(np.isnan(p_rad).all() or np.isnan(p_input).all()):
            (valid_indices,) = np.where((p_input != 0) & ~np.isnan(p_input))
            rad_input_frac[valid_indices] = (
                p_rad[valid_indices] / p_input[valid_indices]
            )
        if ~(np.isnan(p_rad).all() or np.isnan(p_loss).all()):
            (valid_indices,) = np.where((p_loss != 0) & ~np.isnan(p_loss))
            rad_loss_frac[valid_indices] = p_rad[valid_indices] / p_loss[valid_indices]

        output = {
            "p_ecrh": p_ecrh,
            "p_lh": p_lh,
            "p_icrf": p_icrf,
            "p_nbi": p_nbi,
            "rad_input_frac": rad_input_frac,
            "rad_loss_frac": rad_loss_frac,
            "p_input": p_input,
        }
        return output

    @staticmethod
    @physics_method(columns=["p_oh"], tokamak=Tokamak.EAST)
    def get_p_ohm(params: PhysicsMethodParams):
        r"""
        Calculate the ohmic heating power from the loop voltage, inductive voltage, and
        plasma current.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the ohmic heating power (`p_oh`).

        References
        -------
        - original source: [get_P_ohm.m](https://github.com/MIT-PSFC/disruption-py
        /blob/matlab/EAST/get_P_ohm.m)
        """
        # Get raw signals
        try:
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
        except mdsExceptions.MdsException:
            params.logger.warning("Failed to get necessary signals to compute p_ohm")
            return {"p_oh": [np.nan]}
        sign_ip = np.sign(sum(ip))
        dipdt = np.gradient(ip, ip_time)
        # TODO: switch to causal boxcar smoothing
        dipdt_smoothed = matlab_smooth(dipdt, 11)  # Use 11-point boxcar smoothing

        # Interpolate fetched signals to the requested timebase
        vloop = interp1(vloop_time, vloop, params.times, kind="linear", bounds_error=0)
        li = interp1(li_time, li, params.times, kind="linear", bounds_error=0)
        ip = interp1(ip_time, ip, params.times, kind="linear", bounds_error=0)
        dipdt_smoothed = interp1(
            ip_time, dipdt_smoothed, params.times, kind="linear", bounds_error=0
        )

        # Calculate p_ohm
        # TODO: replace 1.85 with fetched r0 signal
        inductance = 4e-7 * np.pi * 1.85 * li / 2  # For EAST use R0 = 1.85 m
        v_inductive = -inductance * dipdt_smoothed
        v_resistive = vloop - v_inductive
        (v_resistive_indices,) = np.where(v_resistive * sign_ip < 0)
        v_resistive[v_resistive_indices] = 0
        p_ohm = ip * v_resistive

        return {"p_oh": p_ohm}

    @staticmethod
    @physics_method(
        columns=[
            "n_equal_1_mode",
            "n_equal_1_phase",
            "n_equal_1_normalized",
            "rmp_n_equal_1",
            "rmp_n_equal_1_phase",
        ],
        tokamak=Tokamak.EAST,
    )
    def get_n_equal_1_data(params: PhysicsMethodParams):
        """
        Compute the amplitude and phase of the n=1 Fourier component of the net
        saddle signals (total saddle signals minus the calculated pickup from
        the RMP coils).

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `n_equal_1_mode`: Amplitude of n=1 Fourier component of saddle signals
            after subtracting the calculated pickup from the RMP coil currents [T].
            - `n_equal_1_phase`: Toroidal phase of the n=1 mode [rad].
            - `n_equal_1_normalized`: `n_equal_1_mode` normalized to the toroidal magnetic
            field (`btor`).
            - `rmp_n_equal_1`: Amplitude of the n=1 Fourier component of the calculated
            pickup of the RMP coils on the saddle signals [T].
            - `rmp_n_equal_1_phase`: toroidal phase of the n=1 rmp pickup [rad].

        References
        -------
        - original source: [get_n_equal_1_data.m](https://github.com/MIT-PSFC/disruption-
        py/blob/matlab/EAST/get_n_equal_1_data.m), [get_rmp_and_saddle_signals.m](https://
        github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_rmp_and_saddle_signals.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411)
        """
        n_equal_1_mode = [np.nan]
        n_equal_1_phase = [np.nan]
        n_equal_1_normalized = [np.nan]
        rmp_n_equal_1 = [np.nan]
        rmp_n_equal_1_phase = [np.nan]

        # Get the rmp coil currents
        # Translated from get_rmp_and_saddle_signals.m
        (rmptime,) = params.mds_conn.get_dims(r"\irmpu1", tree_name="east")
        rmp = np.full((len(rmptime), 16), np.nan)
        for i in range(8):
            # Get irmpu1 to irmpu8
            signal = params.mds_conn.get_data(rf"\irmpu{i+1}", tree_name="east")
            if len(signal) == len(rmptime):
                rmp[:, i] = signal
            # Get irmpl1 to irmpl8
            signal = params.mds_conn.get_data(rf"\irmpl{i+1}", tree_name="east")
            if len(signal) == len(rmptime):
                rmp[:, i + 8] = signal
        # Get saddle coil signals
        (saddletime,) = params.mds_conn.get_dims(r"\sad_pa", tree_name="east")
        saddle = np.full((len(saddletime), 8), np.nan)
        saddle_nodes = [
            r"\sad_pa",
            r"\sad_bc",
            r"\sad_de",
            r"\sad_fg",
            r"\sad_hi",
            r"\sad_jk",
            r"\sad_lm",
            # r"\sad_no",   # not operational in 2015
        ]
        for i, node in enumerate(saddle_nodes[:7]):
            try:
                saddle[:, i] = params.mds_conn.get_data(node, tree_name="east")
            except mdsExceptions.MdsException:
                saddle[:, i] = 0
        sad_lo = params.mds_conn.get_data(r"\sad_lo", tree_name="east")
        sad_lm = params.mds_conn.get_data(r"\sad_lm", tree_name="east")
        saddle[:, 7] = sad_lo - sad_lm

        # Calculate RMP n=1 Fourier component amplitude and phase (on the timebase
        # of the saddle signals)
        coeff_matrix = EastUtilMethods.load_rmp_saddle_coeff_matrix()
        rmp_pickup = np.transpose(np.matmul(coeff_matrix, np.transpose(rmp)))
        # Interpolate rmp_pickup onto the saddle time timebase
        rmp_pickup = interp1(rmptime, rmp_pickup, saddletime, axis=0)

        # Calculate fast Fourier transforms of the RMP-induced signals, and get
        # mode amplitudes and phases
        # Take FFT along 2nd dimension (phi)
        rmp_fft_output = scipy.fft.fft(rmp_pickup, axis=1)
        amplitude = abs(rmp_fft_output) / len(rmp_fft_output[0])
        amplitude[:, 1:] *= 2  # TODO: Why?
        phase = np.arctan2(np.imag(rmp_fft_output), np.real(rmp_fft_output))
        # Only want n=1 Fourier component
        rmp_n_equal_1 = amplitude[:, 1]
        rmp_n_equal_1_phase = phase[:, 1]
        # Interpolate onto the requested timebase
        # TODO: figure out how to interpolate phase
        rmp_n_equal_1 = interp1(saddletime, rmp_n_equal_1, params.times)
        rmp_n_equal_1_phase = interp1(saddletime, rmp_n_equal_1_phase, params.times)

        # The saddle signals can include direct pickup from the RMP coils, when the
        # RMP coils are active.  This direct pickup must be subtracted from the
        # saddle signals to leave just the plasma response and baseline drifts.
        saddle -= rmp_pickup
        # Calculate fast Fourier transforms and get mode amplitudes and phases
        # Take FFT along 2nd dimension (phi)
        saddle_fft_output = scipy.fft.fft(saddle, axis=1)
        amplitude = abs(saddle_fft_output) / len(saddle_fft_output[0])
        amplitude[:, 1:] *= 2  # TODO: Why?
        phase = np.arctan2(np.imag(saddle_fft_output), np.real(saddle_fft_output))
        # Only want n=1 Fourier component
        n_equal_1_mode = amplitude[:, 1]
        n_equal_1_phase = phase[:, 1]
        # Interpolate onto the requested timebase
        n_equal_1_mode = interp1(saddletime, n_equal_1_mode, params.times)
        n_equal_1_phase = interp1(saddletime, n_equal_1_phase, params.times)

        # Next, get the toroidal magnetic field and compute n_equal_1_normalized
        btor = EastPhysicsMethods.get_btor(params)["btor"]
        n_equal_1_normalized = n_equal_1_mode / abs(btor)

        return {
            "n_equal_1_mode": n_equal_1_mode,
            "n_equal_1_phase": n_equal_1_phase,
            "n_equal_1_normalized": n_equal_1_normalized,
            "rmp_n_equal_1": rmp_n_equal_1,
            "rmp_n_equal_1_phase": rmp_n_equal_1_phase,
        }

    @staticmethod
    @physics_method(columns=["btor"], tokamak=Tokamak.EAST)
    def get_btor(params: PhysicsMethodParams):
        r"""
        Calculate the toroidal magnetic field signal.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the toroidal magnetic field signal (`btor`).

        References
        -------
        - original sources: [get_n_equal_1_data.m](https://github.com/MIT-PSFC/disr
        uption-py/blob/matlab/EAST/get_n_equal_1_data.m), [get_mirnov_std.m](https:
        //github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_mirnov_std.m),
        [get_n1rms_n2rms.m](https://github.com/MIT-PSFC/disruption-py/blob/matlab/
        EAST/get_n1rms_n2rms.m)
        """
        if 60000 <= params.shot_id <= 65165:
            tree = "pcs_east"
        else:
            tree = "eng_tree"
        itf, btor_time = params.mds_conn.get_data_with_dims(r"\it", tree_name=tree)
        btor = (
            (4 * np.pi * 1e-7) * itf * (16 * 130) / (2 * np.pi * 1.8)
        )  # about 4,327 amps/tesla
        if params.shot_id < 60000:
            # Btor is constant in time (superconducting magnet).
            # Construct 2-point signal from scalar value
            btor = np.mean(btor)
            btor = [btor, btor]
            btor_time = [0, 1000]

        # Interpolate onto the requested timebase
        btor = interp1(btor_time, btor, params.times)

        return {"btor": btor}

    @staticmethod
    @physics_method(columns=["kappa_area"], tokamak=Tokamak.EAST)
    def get_kappa_area(params: PhysicsMethodParams):
        r"""
        Calculate the plasma's ellipticity (kappa, also known as
        the elongation) using its area and minor radius. It is defined as:

        $$
        \kappa_{area} = \frac{A}{\pi a^2}
        $$

        where $A$ is the plasma cross-sectional area and $a$ is the minor radius.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the `kappa_area` signal.

        References
        -------
        - original source: [get_kappa_area.m](https://github.com/MIT-PSFC/disrup
        tion-py/blob/matlab/EAST/get_kappa_area.m)
        """
        # Get area and aminor from EastEfitMethods
        efit_params = EastEfitMethods.get_efit_parameters(params=params)
        area = efit_params["area"]
        aminor = efit_params["aminor"]
        # Compute kappa_area
        with np.errstate(divide="ignore", invalid="ignore"):
            kappa_area = area / (np.pi * aminor**2)

        return {"kappa_area": kappa_area}

    @staticmethod
    @physics_method(columns=["pkappa_area"], tokamak=Tokamak.EAST)
    def get_pkappa_area(params: PhysicsMethodParams):
        r"""
        Calculate the plasma's ellipticity (kappa) using the area and minor
        radius data from the P-EFIT tree. `kappa_area` is defined as:

        $$
        \kappa_{area} = \frac{A}{\pi a^2}
        $$

        where $A$ is the plasma cross-sectional area and $a$ is the minor radius.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the `pkappa_area` signal.

        References
        -------
        - original source: [get_PEFIT_parameters.m](https://github.com/MIT-PSFC/
        disruption-py/blob/matlab/EAST/get_PEFIT_parameters.m)
        """
        # Get area and aminor from EastEfitMethods
        pefit_params = EastEfitMethods.get_pefit_parameters(params=params)
        area = pefit_params["parea"]
        aminor = pefit_params["paminor"]
        # Compute kappa_area
        with np.errstate(divide="ignore", invalid="ignore"):
            kappa_area = area / (np.pi * aminor**2)

        return {"pkappa_area": kappa_area}

    @staticmethod
    @physics_method(
        columns=["ip_error_rt", "q95_rt", "beta_p_rt", "li_rt", "wmhd_rt"],
        tokamak=Tokamak.EAST,
    )
    def get_pcs_parameters(params: PhysicsMethodParams):
        """
        Retrieve some real-time diagnostic signals that are used in the EAST
        plasma control system (PCS).

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `ip_error_rt`: error between actual and pre-programmed plasma currents.
            - `q95_rt`: q95 calculated by RT-EFIT for PCS (P-EFIT since 2018).
            - `beta_p_rt`: beta_p calculated by RT-EFIT for PCS (P-EFIT since 2018).
            - `li_rt`: li calculated by RT-EFIT for PCS (P-EFIT since 2018).
            - `wmhd_rt`:  Wmhd calculated by RT-EFIT for PCS (P-EFIT since 2018).

        References
        -------
        - original source: [get_pcs.m](https://github.com/MIT-PSFC/disruption-py/
        blob/matlab/EAST/get_pcs.m)
        """
        output = {}

        # Get ip_error_rt, beta_p_rt, li_rt, and wmhd_rt
        signals = {
            "ip_error_rt": r"\lmeip",
            "beta_p_rt": r"\pfsbetap",
            "li_rt": r"\pfsli",
            "wmhd_rt": r"\pfswmhd",
        }
        for name, node in signals.items():
            try:
                signal, timearray = params.mds_conn.get_data_with_dims(
                    node, tree_name="pcs_east"
                )
                signal = interp1(
                    timearray, signal, params.times, kind="linear", bounds_error=0
                )
                output[name] = signal
            except mdsExceptions.MdsException:
                output[name] = [np.nan]

        # Get q95_rt
        try:
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
        except mdsExceptions.MdsException:
            output["q95_rt"] = [np.nan]

        return output

    @staticmethod
    @physics_method(columns=["p_rad_rt", "p_lh_rt", "p_nbi_rt"], tokamak=Tokamak.EAST)
    def get_pcs_power(params: PhysicsMethodParams):
        """
        Calculate some real-time auxiliary heating power signals that are available
        in the EAST plasma control system (PCS).

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `p_rad_rt`: radiated power [W].
            - `p_lh_RT`: lower hybrid power [W].
            - `p_nbi_RT`: neutral beam injection power [W].

        References
        -------
        - original source: [get_pcs.m](https://github.com/MIT-PSFC/disruption-py
        /blob/matlab/EAST/get_pcs.m)
        """
        # Get p_rad_rt
        try:
            p_rad_rt, timearray = params.mds_conn.get_data_with_dims(
                r"\pcprad", tree_name="pcs_east"
            )
            p_rad_rt = interp1(
                timearray, p_rad_rt, params.times, kind="linear", bounds_error=0
            )
        except mdsExceptions.MdsException:
            p_rad_rt = [np.nan]

        # TODO: Verify the outputs, then modify get_heating_power() to make it compatible with this
        # Get p_nbi_rt
        nbi_nodes = {
            "i_nbil_rt": r"\pcnbi1li",
            "v_nbil_rt": r"\pcnbi1lv",
            "i_nbir_rt": r"\pcnbi1ri",
            "v_nbir_rt": r"\pcnbi1rv",
        }
        try:
            nbi_signals = {}
            for name, node in nbi_nodes.items():
                signal, timearray = params.mds_conn.get_data_with_dims(
                    node, tree_name="pefitrt_east"
                )
                signal = interp1(
                    timearray, signal, params.times, kind="linear", bounds_error=0
                )
                nbi_signals[name] = signal
            p_nbi_rt = (
                nbi_signals["i_nbil_rt"] * nbi_signals["v_nbil_rt"]
                + nbi_signals["i_nbir_rt"] * nbi_signals["v_nbir_rt"]
            )
        except mdsExceptions.MdsException:
            p_nbi_rt = [np.nan]

        # Get p_lh_rt
        lh_nodes = {
            "p_lh_46_inj_rt": r"\pcplhi",
            "p_lh_46_ref_rt": r"\pcplhr",
            "p_lh_245_inj_rt": r"\pcplhi2",
            "p_lh_245_ref_rt": r"\pcplhr2",
        }
        try:
            lh_signals = {}
            for name, node in lh_nodes.items():
                signal, timearray = params.mds_conn.get_data_with_dims(
                    node, tree_name="pefitrt_east"
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
        except mdsExceptions.MdsException:
            p_lh_rt = [np.nan]

        # Q.P. Yuan:
        # There is no signals from ICRF or ECRH connected to PCS. And I checked the
        # signals for NBI, there is no effective value, just noise. We will make
        # the NBI signals available in this EAST campaign. The PLHI2 and PLHR2 for
        # 2.45G LHW system have signals but are not calibrated yet.

        # Get p_ohm
        # "Note, I'm not bothering to calculate a real time version" -- Whoever who wrote this
        # p_ohm = EastPhysicsMethods.get_p_ohm(params)['p_ohm']

        # Compute total power and radiation fractions
        # p_input_rt = p_ohm + p_lh_rt + p_icrf_rt + p_ecrh_rt + p_nbi_rt
        # rad_input_frac_rt = p_rad_rt / p_input_rt
        # rad_loss_frac_rt = p_rad_rt / (p_input_rt - dwmhd_dt)

        return {
            "p_rad_rt": p_rad_rt,
            "p_lh_rt": p_lh_rt,
            "p_nbi_rt": p_nbi_rt,
        }

    @staticmethod
    @physics_method(columns=["prad_peaking"], tokamak=Tokamak.EAST)
    def get_prad_peaking(params: PhysicsMethodParams):
        """
        Calculates the peaking factor of the profiles of radiated
        power measured by the AXUV arrays on EAST. It is defined as the
        ratio of the average of the centralmost 6 channels to the average
        of all of the non-divertor-viewing channels (core-to-average).

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the peaking factor for radiated power (`prad_peaking`).

        References
        -------
        - original source: [get_prad_peaking.m](https://github.com/MIT-PSFC/disruption-py/blob
        /matlab/EAST/get_prad_peaking.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411), #[451](https:
        //github.com/MIT-PSFC/disruption-py/pull/451)
        """
        xuv, xuvtime = EastPhysicsMethods._get_raw_axuv_data(params)

        # Get calibration factors
        calib_factors = EastUtilMethods.get_axuv_calib_factors()
        # Apply calibration factors
        for i in range(64):
            # Original script has Del_r(ichan) which should be a mistake
            xuv[:, i] *= (
                calib_factors["Fac1"][i % 16]
                * calib_factors["Fac2"][i // 16][i % 16]
                * calib_factors["Fac3"][i // 16]
                * calib_factors["Fac4"]
                * 2
                * np.pi
                * calib_factors["Maj_R"]
                * calib_factors["Del_r"][i]
                # * calib_factors["Fac5"]
            )
        # Correction for bad channels (from Duan Yanmin's program)
        # get_radiated_power does these 2 lines first and then apply geometric correction.
        xuv[:, 11] = 0.5 * (xuv[:, 10] + xuv[:, 12])
        # xuv[:, 35] = 0.5 * (xuv[:, 34] + xuv[:, 35])

        # Define the core chords to be #28 to #37 (centermost 10 chords), and
        # define the non-divertor chords to be #09 to #56.  (Chords #1-8 view the
        # lower divertor region, and chords #57-64 view the upper divertor region.)
        ch_core = np.zeros(64)
        ch_core[29:35] = 1 / 6
        ch_all = np.zeros(64)
        ch_all[8:56] = 1 / 48

        prad_peaking = np.full((len(xuvtime)), np.nan)
        # TODO: find better way to do this for loop
        for itime in range(len(xuvtime)):
            prad_core = np.dot(xuv[itime, :], ch_core)
            prad_all = np.dot(xuv[itime, :], ch_all)
            prad_peaking[itime] = prad_core / prad_all

        # Interpret to the requested timebase
        prad_peaking = interp1(xuvtime, prad_peaking, params.times)

        return {"prad_peaking": prad_peaking}

    @staticmethod
    @physics_method(
        columns=["mirnov_std", "mirnov_std_normalized"], tokamak=Tokamak.EAST
    )
    def get_mirnov_std(params: PhysicsMethodParams):
        """
        Fetch the signal from a single Mirnov sensor, then compute the rolling
        standard deviation.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `mirnov_std`: rolling stdev of a mirnov sensor.
            - `mirnov_std_normalized`: `mirnov_std` normalized to the toroidal
            magnetic field (`btor`).

        References
        -------
        - original source: [get_mirnov_std.m](https://github.com/MIT-PSFC/disrupt
        ion-py/blob/matlab/EAST/get_mirnov_std.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411)
        """
        mirnov_std = np.full(len(params.times), np.nan)
        mirnov_std_normalized = [np.nan]

        # Get the toroidal magnetic field.
        btor = EastPhysicsMethods.get_btor(params)["btor"]

        # Get the Mirnov signal from \cmp1t (5 MHz)
        time_window = 0.001
        bp_dot, bp_dot_time = params.mds_conn.get_data_with_dims(
            r"\cmp1t", tree_name="east"
        )  # [T/s], [s]
        for i, time in enumerate(params.times):
            (indices,) = np.where(
                ((time - time_window) < bp_dot_time) & (bp_dot_time < time)
            )
            mirnov_std[i] = np.nanstd(bp_dot[indices], ddof=1)

        mirnov_std_normalized = mirnov_std / abs(btor)

        return {
            "mirnov_std": mirnov_std,
            "mirnov_std_normalized": mirnov_std_normalized,
        }

    @staticmethod
    @physics_method(
        columns=["n1rms", "n2rms", "n1rms_normalized", "n2rms_normalized"],
        tokamak=Tokamak.EAST,
    )
    def get_n1rms_n2rms(params: PhysicsMethodParams):
        """
        Retrieve data of the Mirnov sensor array and compute the n=1 and 2 Fourier
        mode amplitudes. Then calculate the rolling standard deviation of the amplitudes
        and normalize them to the toroidal magnetic field (`btor`).

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the MDS connection and shot information.

        Returns
        -------
        dict
            A dictionary containing the following keys:

            - `n1rms`: rolling standard deviation of the n=1 mode from the Mirnov array [T].
            - `n2rms`: rolling standard deviation the n=2 mode from the Mirnov array [T].
            - `n1rms_normalized`: `n1rms` normalized to `btor` [dimensionless].
            - `n2rms_normalized`: `n2rms` normalized to `btor` [dimensionless].

        References
        -------
        - original source: [get_n1rms_n2rms.m](https://github.com/MIT-PSFC/disruption-py/blob
        /matlab/EAST/get_n1rms_n2rms.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411)
        """
        n1rms = np.full(len(params.times), np.nan)
        n2rms = np.full(len(params.times), np.nan)
        n1rms_normalized = [np.nan]
        n2rms_normalized = [np.nan]

        (mirtime,) = params.mds_conn.get_dims(r"\mitab2", tree_name="east")
        mir = np.full((len(mirtime), 16), np.nan)
        mir_nodes = [
            r"\mitab2",
            r"\mitbc2",
            r"\mitcd2",
            r"\mitde2",
            r"\mitef2",
            r"\mitfg2",
            r"\mitgh2",
            r"\mithi2",
            r"\mitij2",
            r"\mitjk2",
            r"\mitkl2",
            r"\mitlm2",
            r"\mitmn2",
            r"\mitno2",
            r"\mitop2",
            r"\mitpa2",
        ]
        for i, node in enumerate(mir_nodes):
            try:
                mir[:, i] = params.mds_conn.get_data(node, tree_name="east")
            except mdsExceptions.MdsException:
                continue

        # Get btor data
        btor = EastPhysicsMethods.get_btor(params)["btor"]

        # Compute the n=1 and n=2 mode signals from the Mirnov array signals
        mir_fft_output = scipy.fft.fft(mir, axis=1)
        amplitude = abs(mir_fft_output) / len(mir_fft_output[0])
        amplitude[:, 1:] *= 2  # TODO: Why?
        phase = np.arctan2(np.imag(mir_fft_output), np.real(mir_fft_output))
        n1 = amplitude[:, 1] * np.cos(phase[:, 1])
        n2 = amplitude[:, 2] * np.cos(phase[:, 2])

        # Calculate n1rms & n2rms
        time_window = 0.001
        for i, time in enumerate(params.times):
            (indices,) = np.where(
                (time - time_window < mirtime) & (mirtime < time + time_window)
            )
            n1rms[i] = np.nanstd(n1[indices], ddof=1)
            n2rms[i] = np.nanstd(n2[indices], ddof=1)

        # Calculate the normalized signals
        n1rms_normalized = n1rms / abs(btor)
        n2rms_normalized = n2rms / abs(btor)

        return {
            "n1rms": n1rms,
            "n2rms": n2rms,
            "n1rms_normalized": n1rms_normalized,
            "n2rms_normalized": n2rms_normalized,
        }

    @staticmethod
    @physics_method(columns=["h98"], tokamak=Tokamak.EAST)
    def get_h98(params: PhysicsMethodParams):
        """
        Get the H98y2 energy confinement time parameter.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the H98y2 energy confinement time (`h98`).

        References
        -------
        - original source: [get_h98.m](https://github.com/MIT-PSFC/disruption-py/blob/
        matlab/EAST/get_h98.m)
        """
        h98_y2 = [np.nan]

        h98_y2, h98_y2_time = params.mds_conn.get_data_with_dims(
            r"\h98_mhd", tree_name="energy_east"
        )

        # Interpolate the signal onto the requested timebase
        h98_y2 = interp1(h98_y2_time, h98_y2, params.times)

        return {"h98": h98_y2}

    @staticmethod
    def _get_efit_gaps(params: PhysicsMethodParams, tree: str = "_efit_tree"):
        """
        Hidden method to calculate the EFIT and P-EFIT gaps
        """
        upper_gap = [np.nan]
        lower_gap = [np.nan]

        # Get plasma boundary data
        data, efittime = params.mds_conn.get_data_with_dims(
            r"\top.results.geqdsk:bdry", tree_name=tree
        )
        # Convert the order of indices to MATLAB order
        # MATLAB (0, 1, 2) -> Python (2, 1, 0)
        data = np.transpose(data, [2, 1, 0])
        xcoords, ycoords = data

        # Get first wall geometry data
        xfirstwall = params.mds_conn.get_data(
            r"\top.results.geqdsk:xlim", tree_name=tree
        )
        yfirstwall = params.mds_conn.get_data(
            r"\top.results.geqdsk:ylim", tree_name=tree
        )
        seed = np.ones((len(xcoords), 1))
        xfirstwall = np.reshape(xfirstwall, (-1, 1))
        yfirstwall = np.reshape(yfirstwall, (-1, 1))
        xfirstwall_mat = np.tile(xfirstwall, (1, len(efittime)))
        yfirstwall_mat = np.tile(yfirstwall, (1, len(efittime)))

        # Calculate upper & lower gaps
        index_upperwall, _ = np.where(yfirstwall > 0.6)
        index_lowerwall, _ = np.where(yfirstwall < -0.6)

        xupperwall = xfirstwall_mat[index_upperwall, :]
        xupperwall_mat = np.reshape(
            np.kron(xupperwall, seed), (len(seed), -1, len(efittime))
        )
        xlowerwall = xfirstwall_mat[index_lowerwall, :]
        xlowerwall_mat = np.reshape(
            np.kron(xlowerwall, seed), (len(seed), -1, len(efittime))
        )

        yupperwall = yfirstwall_mat[index_upperwall, :]
        yupperwall_mat = np.reshape(
            np.kron(yupperwall, seed), (len(seed), -1, len(efittime))
        )
        ylowerwall = yfirstwall_mat[index_lowerwall, :]
        ylowerwall_mat = np.reshape(
            np.kron(ylowerwall, seed), (len(seed), -1, len(efittime))
        )

        xupperplasma_mat = np.reshape(
            np.tile(xcoords, (len(xupperwall), 1)), (-1, len(xupperwall), len(efittime))
        )
        yupperplasma_mat = np.reshape(
            np.tile(ycoords, (len(yupperwall), 1)), (-1, len(yupperwall), len(efittime))
        )
        xlowerplasma_mat = np.reshape(
            np.tile(xcoords, (len(xlowerwall), 1)), (-1, len(xlowerwall), len(efittime))
        )
        ylowerplasma_mat = np.reshape(
            np.tile(ycoords, (len(ylowerwall), 1)), (-1, len(ylowerwall), len(efittime))
        )

        uppergap_mat = np.sqrt(
            np.square(xupperplasma_mat - xupperwall_mat)
            + np.square(yupperplasma_mat - yupperwall_mat)
        )
        lowergap_mat = np.sqrt(
            np.square(xlowerplasma_mat - xlowerwall_mat)
            + np.square(ylowerplasma_mat - ylowerwall_mat)
        )
        upper_gap = uppergap_mat.min(axis=(0, 1))
        lower_gap = lowergap_mat.min(axis=(0, 1))

        # Interpolate to the requested timebase
        upper_gap = interp1(efittime, upper_gap, params.times)
        lower_gap = interp1(efittime, lower_gap, params.times)

        return {"upper_gap": upper_gap, "lower_gap": lower_gap}

    @staticmethod
    @physics_method(
        columns=["upper_gap", "lower_gap", "pupper_gap", "plower_gap"],
        tokamak=Tokamak.EAST,
    )
    def get_efit_gaps(params: PhysicsMethodParams):
        """
        Calculate the upper and lower gaps from the EFIT and P-EFIT
        information on plasma boundary and first wall geometry.

        Parameters
        ----------
        params : PhysicsMethodParams
            Parameters containing MDS connection and shot information

        Returns
        -------
        dict
            A dictionary containing the upper and lower gaps from EFIT (`upper_gap`
            and `lower_gap) and P-EFIT (`pupper_gap` and `plower_gap`) data.

        References
        -------
        - original sources: [get_EFIT_gaps.m](https://github.com/MIT-PSFC/disruption-
        py/blob/matlab/EAST/get_EFIT_gaps.m), [get_PEFIT_gaps.m](https://github.com/M
        IT-PSFC/disruption-py/blob/matlab/EAST/get_PEFIT_gaps.m)
        - pull requests: #[411](https://github.com/MIT-PSFC/disruption-py/pull/411)
        """
        efit_gaps = EastPhysicsMethods._get_efit_gaps(params=params, tree="_efit_tree")
        pefit_gaps = EastPhysicsMethods._get_efit_gaps(params=params, tree="pefit_east")

        return {
            "upper_gap": efit_gaps["upper_gap"],
            "lower_gap": efit_gaps["lower_gap"],
            "pupper_gap": pefit_gaps["upper_gap"],
            "plower_gap": pefit_gaps["lower_gap"],
        }
