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

        Original author: Robert Granetz, Apr 2016
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
        https://github.com/MIT-PSFC/disruption-py/blob/matlab/EAST/get_Ip_parameters.m

        Original author: Robert Granetz, Dec 2015
        Last major update: 11/19/24 by William Wei
        """
