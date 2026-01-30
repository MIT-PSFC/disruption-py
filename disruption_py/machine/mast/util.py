#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np
import xarray as xr

from disruption_py.inout.xarray_backend import XarrayConnection


class MastUtilMethods:
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
    def retrieve_ip(conn: XarrayConnection, shot_id: int):
        """
        Read in the measured plasma current, Ip.

        Parameters
        ----------
        conn : XarrayConnection
            Connection to S3 bucket.
        shot_id : int
            Shot number.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Plasma current [A], time base of plasma current [s].
        """
        file_name = conn.get_shot_file_path(shot_id)
        ds = xr.open_zarr(file_name, group="summary")
        ip = ds["ip"].values
        ip_time = ds["time"].values
        ip = MastUtilMethods.subtract_ip_baseline_offset(ip, ip_time)
        return ip, ip_time

    @staticmethod
    def retrieve_efit_time(conn: XarrayConnection, shot_id: int):
        """
        Read in the EFIT time base.

        Parameters
        ----------
        conn : XarrayConnection
            Connection to S3 bucket.
        shot_id : int
            Shot number.

        Returns
        -------
        np.ndarray
            EFIT time base [s].
        """
        file_name = conn.get_shot_file_path(shot_id)
        ds = xr.open_zarr(file_name, group="equilibrium")
        efit_time = ds["time"].values
        return efit_time
