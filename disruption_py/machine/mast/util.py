#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np

from disruption_py.core.physics_method.errors import (
    CalculationError,
    MismatchCalculationError,
)
from disruption_py.core.utils.math import interp1
from disruption_py.inout.xr import XarrayDataConnection


class MastUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

    @staticmethod
    def get_signal(conn: XarrayDataConnection, shot_id: int, path: str) -> np.ndarray:
        """
        Read a signal, squeezed but never 0-dimensional.

        `XarrayConnection.get_data` returns a length-1 NaN array for a variable that is
        absent from the shot. Squeezing that sentinel collapses it to a 0-d array, which
        then fails deep inside the numerics rather than at the point of the missing data.

        Parameters
        ----------
        conn : XarrayConnection
            Connection to the data store.
        shot_id : int
            Shot number.
        path : str
            Path to the variable in the shot's data tree.

        Returns
        -------
        np.ndarray
            The signal, with singleton dimensions dropped but at least 1-D.
        """
        return np.atleast_1d(conn.get_data(shot_id, path).squeeze())

    @staticmethod
    def require_signal(conn: XarrayConnection, shot_id: int, path: str) -> np.ndarray:
        """
        Read a signal that a calculation cannot proceed without.

        Parameters
        ----------
        conn : XarrayConnection
            Connection to the data store.
        shot_id : int
            Shot number.
        path : str
            Path to the variable in the shot's data tree.

        Returns
        -------
        np.ndarray
            The signal, with singleton dimensions dropped but at least 1-D.

        Raises
        ------
        CalculationError
            If the signal is absent from the shot or holds no finite samples.
        """
        signal = MastUtilMethods.get_signal(conn, shot_id, path)
        if not np.any(np.isfinite(signal)):
            raise CalculationError(f"No valid data for signal: {path}")
        return signal

    @staticmethod
    def require_time_base(
        conn: XarrayConnection, shot_id: int, path: str
    ) -> np.ndarray:
        """
        Read a time base, which needs two samples before anything can be
        interpolated or differentiated against it.

        Parameters
        ----------
        conn : XarrayConnection
            Connection to the data store.
        shot_id : int
            Shot number.
        path : str
            Path to the time variable in the shot's data tree.

        Returns
        -------
        np.ndarray
            The time base.

        Raises
        ------
        CalculationError
            If the time base is absent, all-NaN, or shorter than two samples.
        """
        time = MastUtilMethods.require_signal(conn, shot_id, path)
        if time.ndim != 1 or time.size < 2:
            raise CalculationError(
                f"Time base {path} needs at least 2 samples, got shape {time.shape}"
            )
        return time

    @staticmethod
    def require_aligned(path: str, signal: np.ndarray, time: np.ndarray) -> np.ndarray:
        """
        Check a signal is sampled on the time base it will be interpolated against.

        The last axis is the time axis for both 1-D traces and 2-D (channel, time)
        profiles, so this covers both.

        Parameters
        ----------
        path : str
            Path of the signal, used for the error message.
        signal : np.ndarray
            The signal to check.
        time : np.ndarray
            The time base the signal should be sampled on.

        Returns
        -------
        np.ndarray
            The signal, unchanged.

        Raises
        ------
        CalculationError
            If the signal's time axis does not match the length of the time base.
        """
        if signal.shape[-1] != len(time):
            raise CalculationError(
                f"{path} has {signal.shape[-1]} samples "
                f"but its time base has {len(time)}"
            )
        return signal

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
        # a missing signal squeezes down to a 0-d array, which has no len()
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

        if np.isnan(y).all() or len(x) < 2:
            # if all y are NaN (is a missing signal)
            # or if x has only a single number
            # just return array of NaNs with same shape as x_new
            return np.full_like(x_new, np.nan)

        if len(x) != len(y):
            raise MismatchCalculationError(f"len(x) = {len(x)} vs. len(y) = {len(y)}")

        return interp1(x, y, x_new)
