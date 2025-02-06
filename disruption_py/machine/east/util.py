"""
Module for helper, not physics, methods.
"""

import numpy as np
import scipy

from disruption_py.inout.mds import MDSConnection


class EastUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

    @staticmethod
    def subtract_ip_baseline_offset(ip, ip_time):
        """
        Subtract the baseline offset from the plasma current signal.

        Parameters
        ----------
        ip : np.ndarray
            Plasma current [A].
        ip_time : np.ndarray
            Time base of plasma current [s].

        Returns
        -------
        np.ndarray
            Offset-subtracted plasma current [A].
        """
        # Get indices of time before any PF supplies turn on
        (base_indices,) = np.where(ip_time <= -5.8)
        # Subtract offset
        if len(base_indices) > 0:
            baseline = sum(ip[base_indices]) / len(base_indices)
            ip -= baseline
        return ip

    @staticmethod
    def retrieve_ip(mds_conn: MDSConnection, shot_id: int):
        """
        Read in the measured plasma current, Ip. There are several different
        measurements of Ip: IPE, IPG, IPM (all in the EAST tree), and PCRL01
        (in the PCS_EAST tree). At various times in the history of EAST, there
        have been problems with all of these measurements, such as broken sensors,
        inverted signals, and shifted timebases. I think the most reliable one is
        PCRL01, which is the one used by the Plasma Control System (PCS) for feedback
        control. So that is the one I will use for the disruption warning database.

        Parameters
        ----------
        mds_conn : MDSConnection
            Connection to MDSplus server.
        shot_id : int
            Shot number.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Plasma current [A], time base of plasma current [s].
        """

        ip, ip_time = mds_conn.get_data_with_dims(
            r"\pcrl01", tree_name="pcs_east"
        )  # [A], [s]

        # For shots before year 2014, the PCRL01 timebase needs to be shifted
        # by 17.0 ms
        if shot_id < 44432:
            ip_time -= 0.0170

        # High-frequency noise spikes on some shots can cause a problem with the
        # time derivative and other computations.  Use a median filter to reduce
        # the problem.
        ip = scipy.signal.medfilt(ip, 5)  # Remove noise spikes with median filter

        # Subtract baseline offset
        ip = EastUtilMethods.subtract_ip_baseline_offset(ip, ip_time)
        return ip, ip_time
