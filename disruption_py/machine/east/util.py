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

    @staticmethod
    def get_axuv_calib_factors() -> dict:
        """
        Return the calibration factors for the AXUV arrays

        Used by get_radiated_power and get_prad_peaking.

        - fac_1: (unknown)
        - fac_2: factors of Amp.Gain
        - fac_3: cross calibration factors between arrays
        - fac_4: unit convert
        - fac_5: corrected factor by cross calibration with foil bolometer
        - maj_r: major radius (of the machine cross-section?)
        - del_r: (unknown)
        """
        fac_1 = (
            np.array(
                [
                    1.3681,
                    1.3429,
                    1.3215,
                    1.3039,
                    1.2898,
                    1.2793,
                    1.2723,
                    1.2689,
                    1.2689,
                    1.2723,
                    1.2793,
                    1.2898,
                    1.3039,
                    1.3215,
                    1.3429,
                    1.3681,
                ]
            )
            * 1e4
        )
        fac_2 = np.array(
            [
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            ]
        )
        fac_3 = np.array([1, 1, 1, 1])
        del_r = (
            np.array(
                [
                    3.6,
                    3.6,
                    3.5,
                    3.4,
                    3.4,
                    3.3,
                    3.3,
                    3.2,
                    3.2,
                    3.1,
                    3.1,
                    3.0,
                    3.0,
                    2.9,
                    2.9,
                    2.8,
                    2.9,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.6,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.7,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.8,
                    2.9,
                    2.8,
                    2.9,
                    2.9,
                    3.0,
                    3.0,
                    3.1,
                    3.1,
                    3.2,
                    3.2,
                    3.3,
                    3.3,
                    3.4,
                    3.4,
                    3.5,
                    3.6,
                    3.6,
                ]
            )
            * 0.01
        )
        # Use original names in the MATLAB script
        return {
            "Fac1": fac_1,
            "Fac2": fac_2,
            "Fac3": fac_3,
            "Fac4": 1e-3,
            "Fac5": 2.5,
            "Maj_R": 1.85,
            "Del_r": del_r,
        }
