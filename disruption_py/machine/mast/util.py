#!/usr/bin/env python3

"""
Module for helper, not physics, methods.
"""

import numpy as np

from disruption_py.core.physics_method.errors import (
    CalculationError,
    MismatchCalculationError,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1


class MastUtilMethods:
    """
    A class of helper methods that might fetch and compute data from MDSplus
    but are not physics methods.
    """

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
    def retrieve_ip(params: PhysicsMethodParams):
        """
        Read in the measured plasma current, Ip.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection and shot id.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Plasma current [A], time base of plasma current [s].
        """
        ip = params.get_data("summary/ip")
        ip_time = params.get_data("summary/time")
        return ip, ip_time

    @staticmethod
    def retrieve_efit_time(params: PhysicsMethodParams):
        """
        Read in the EFIT time base.

        Parameters
        ----------
        params : PhysicsMethodParams
            Per-shot Xarray data connection.

        Returns
        -------
        np.ndarray
            EFIT time base [s].
        """
        efit_time = params.get_data("equilibrium/time")
        return efit_time

    @staticmethod
    def thomson_rho(params: PhysicsMethodParams, r_ts: np.ndarray):
        """
        Normalised effective radius rho = |R_TS - R_mag| / a_minor for each Thomson
        scattering channel.

        The Thomson channels sit at fixed major radius, so time-averaged equilibrium
        values are used. The averages are returned alongside rho because callers also
        need them to bound profile fits in machine coordinates.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection and shot id.
        r_ts : np.ndarray
            Major radius of each Thomson scattering channel [m].

        Returns
        -------
        tuple[np.ndarray, float, float]
            Normalised effective radius per channel, time-averaged major radius of
            the magnetic axis [m], time-averaged minor radius [m].

        Raises
        ------
        CalculationError
            If the equilibrium magnetic axis or minor radius is unavailable.
        """
        r_mag = params.get_data("equilibrium/magnetic_axis_r")
        a_minor = params.get_data("equilibrium/minor_radius")
        r_mag_mean = np.nanmean(r_mag)
        a_minor_mean = np.nanmean(a_minor)
        if (
            not np.isfinite(r_mag_mean)
            or not np.isfinite(a_minor_mean)
            or a_minor_mean <= 0
        ):
            raise CalculationError(
                "Cannot compute rho for Thomson scattering channels: "
                "equilibrium magnetic_axis_r or minor_radius unavailable."
            )
        return np.abs(r_ts - r_mag_mean) / a_minor_mean, r_mag_mean, a_minor_mean

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
        # a 0-d array has no len()
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
