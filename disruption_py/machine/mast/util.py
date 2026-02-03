#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np

from disruption_py.inout.xr import XarrayConnection
from disruption_py.core.utils.math import interp1


class MastUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

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
        ip = conn.get_data(shot_id, "summary/ip")  # ensure data is cached
        ip_time = conn.get_data(shot_id, "summary/time")
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
        efit_time = conn.get_data(shot_id, "equilibrium/time")
        return efit_time

    @staticmethod
    def interpolate_1d(x, y, x_new):
        """Safely interpolate 1D data with handling for all-NaN y values.

        Parameters
        ----------
        x : array_like
            Original x-coordinates of the data points.
        y : array_like
            Original y-coordinates of the data points.
        x_new : array_like
            New x-coordinates where interpolation is desired.

        Returns
        -------
        array_like
            Interpolated y-coordinates corresponding to x_new.
        """
        if len(x) != len(y) and np.isnan(y).all():
            # if all y are NaN (is a missing signal)
            # just return array of NaNs with same shape as x_new
            return np.full_like(x_new, np.nan)
        return interp1(x, y, x_new)
