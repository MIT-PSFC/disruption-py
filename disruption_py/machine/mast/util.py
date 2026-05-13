#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np

from disruption_py.core.utils.math import interp1
from disruption_py.inout.xr import XarrayDataConnection


class MastUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

    @staticmethod
    def retrieve_ip(conn: XarrayDataConnection):
        """
        Read in the measured plasma current, Ip.

        Parameters
        ----------
        conn : XarrayDataConnection
            Per-shot Xarray data connection.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Plasma current [A], time base of plasma current [s].
        """
        ip = conn.get_data("summary/ip")
        ip_time = conn.get_data("summary/time")
        return ip, ip_time

    @staticmethod
    def retrieve_efit_time(conn: XarrayDataConnection):
        """
        Read in the EFIT time base.

        Parameters
        ----------
        conn : XarrayDataConnection
            Per-shot Xarray data connection.

        Returns
        -------
        np.ndarray
            EFIT time base [s].
        """
        efit_time = conn.get_data("equilibrium/time")
        return efit_time

    @staticmethod
    def interpolate_1d(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
        """Safely interpolate 1D data with handling for all-NaN y values.

        Parameters
        ----------
        x : np.ndarray
            Original x-coordinates of the data points.
        y : np.ndarray
            Original y-coordinates of the data points.
        x_new : np.ndarray
            New x-coordinates where interpolation is desired.

        Returns
        -------
        np.ndarray
            Interpolated y-coordinates corresponding to x_new.
        """
        if np.isnan(y).all() or len(x) < 2:
            # if all y are NaN (is a missing signal)
            # or if x has only a single number
            # just return array of NaNs with same shape as x_new
            return np.full_like(x_new, np.nan)

        if len(x) != len(y):
            raise ValueError("x and y must have the same length for interpolation.")

        return interp1(x, y, x_new)
