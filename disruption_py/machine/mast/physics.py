#!/usr/bin/env python3
# pylint: disable=duplicate-code

"""
Physics methods for MAST.
"""

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.inout.xarray import XarrayConnection
from disruption_py.machine.tokamak import Tokamak


class MastPhysicsMethods:
    """
    This class provides methods to retrieve and calculate physics-related data
    for MAST.
    """

    @staticmethod
    @physics_method(
        columns=["ip"],
        tokamak=Tokamak.MAST,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """Get Ip parameters

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing plasma current (`ip`), its time derivative (`dip_dt`),
            programmed plasma current (`ip_prog`), and its time derivative (`dipprog_dt`).
        """
        conn: XarrayConnection = params.mds_conn
        ip = conn.get_data(params.shot_id, "summary/ip")
        ip_prog = conn.get_data(params.shot_id, "pulse_schedule/i_plasma")
        ip_prog_time = conn.get_data(params.shot_id, "pulse_schedule/time")
        magtime = conn.get_data(params.shot_id, "summary/time")

        dip_dt = np.gradient(ip, magtime)
        dipprog_dt = np.gradient(ip_prog, ip_prog_time)

        times = params.times

        ip = MastPhysicsMethods.interpolate_1d(magtime, ip, times)
        ip_prog = MastPhysicsMethods.interpolate_1d(ip_prog_time, ip_prog, times)
        dip_dt = MastPhysicsMethods.interpolate_1d(magtime, dip_dt, times)
        dipprog_dt = MastPhysicsMethods.interpolate_1d(ip_prog_time, dipprog_dt, times)

        return {
            "ip": ip,
            "dip_dt": dip_dt,
            "ip_prog": ip_prog,
            "dipprog_dt": dipprog_dt,
        }

    @staticmethod
    @physics_method(
        columns=["power_nbi", "power_radiated"],
        tokamak=Tokamak.MAST,
    )
    def get_power(params: PhysicsMethodParams):
        """Get power parameters

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing neutral beam injection power (`power_nbi`) and
            radiated power (`power_radiated`).
        """
        conn: XarrayConnection = params.mds_conn

        power_nbi = conn.get_data(params.shot_id, "summary/power_nbi")
        power_radiated = conn.get_data(params.shot_id, "summary/power_radiated")
        base_time = conn.get_data(params.shot_id, "summary/time")

        times = params.times
        power_nbi = MastPhysicsMethods.interpolate_1d(base_time, power_nbi, times)
        power_radiated = MastPhysicsMethods.interpolate_1d(
            base_time, power_radiated, times
        )
        return {"power_nbi": power_nbi, "power_radiated": power_radiated}

    @staticmethod
    @physics_method(
        columns=["total_injected", "inboard_total", "outboard_total"],
        tokamak=Tokamak.MAST,
    )
    def get_gas(params: PhysicsMethodParams):
        """Get gas injection parameters

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing total injected gas (`total_injected`),
            inboard total gas (`inboard_total`), and outboard total gas (`outboard_total`).
        """
        conn: XarrayConnection = params.mds_conn

        total_injected = conn.get_data(params.shot_id, "gas_injection/total_injected")
        inboard_total = conn.get_data(params.shot_id, "gas_injection/inboard_total")
        outboard_total = conn.get_data(params.shot_id, "gas_injection/outboard_total")
        base_time = conn.get_data(params.shot_id, "gas_injection/time")

        times = params.times
        total_injected = MastPhysicsMethods.interpolate_1d(
            base_time, total_injected, times
        )
        inboard_total = MastPhysicsMethods.interpolate_1d(
            base_time, inboard_total, times
        )
        outboard_total = MastPhysicsMethods.interpolate_1d(
            base_time, outboard_total, times
        )
        return {
            "total_injected": total_injected,
            "inboard_total": inboard_total,
            "outboard_total": outboard_total,
        }

    @staticmethod
    @physics_method(
        columns=["t_e_core", "n_e_core"],
        tokamak=Tokamak.MAST,
    )
    def get_ts_parameters(params: PhysicsMethodParams):
        """Get Thomson parameters

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing core electron temperature (`t_e_core`) and
            core electron density (`n_e_core`).
        """
        times = params.times
        conn: XarrayConnection = params.mds_conn

        t_e_core = conn.get_data(params.shot_id, "thomson_scattering/t_e_core")
        n_e_core = conn.get_data(params.shot_id, "thomson_scattering/n_e_core")
        base_time = conn.get_data(params.shot_id, "thomson_scattering/time")

        t_e_core = MastPhysicsMethods.interpolate_1d(base_time, t_e_core, times)
        n_e_core = MastPhysicsMethods.interpolate_1d(base_time, n_e_core, times)
        return {"t_e_core": t_e_core, "n_e_core": n_e_core}

    @staticmethod
    @physics_method(
        columns=["n_e", "dn_dt", "greenwald_fraction"],
        tokamak=Tokamak.MAST,
    )
    def get_densities(params: PhysicsMethodParams):
        r"""
        Calculate electron density, its time derivative, and the Greenwald fraction.

        The Greenwald fraction is the ratio of the measured electron density $n_e$ and
        the Greenwald density limit $n_G$ defined as [^1]:

        $$
        n_G = \frac{I_p}{\pi a^2}
        $$

        where $n_G$ is given in $10^{20} m^{-3}$ and $I_p$ is in MA.

        [^1]: https://wiki.fusion.ciemat.es/wiki/Greenwald_limit

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing electron density (`n_e`), its gradient (`dn_dt`),
            and the Greenwald fraction (`greenwald_fraction`).

        """

        conn: XarrayConnection = params.mds_conn
        n_e = conn.get_data(params.shot_id, "summary/line_average_n_e")
        t_n = conn.get_data(params.shot_id, "summary/time")
        ip = conn.get_data(params.shot_id, "summary/ip")
        t_ip = conn.get_data(params.shot_id, "summary/time")
        a_minor = conn.get_data(params.shot_id, "equilibrium/minor_radius")
        t_a = conn.get_data(params.shot_id, "equilibrium/time")

        return MastPhysicsMethods._get_densities(
            params.times, n_e, t_n, ip, t_ip, a_minor, t_a
        )

    @staticmethod
    def _get_densities(times, n_e, t_n, ip, t_ip, a_minor, t_a):
        """
        Calculate electron density, its time derivative, and the Greenwald fraction.

        Parameters
        ----------
        times : array_like
            Time points at which to interpolate the densities.
        n_e : array_like
            Electron density values.
        t_n : array_like
            Corresponding time values for electron density.
        ip : array_like
            Plasma current values.
        t_ip : array_like
            Corresponding time values for plasma current.
        a_minor : array_like
            Minor radius values.
        t_a : array_like
            Corresponding time values for minor radius.

        Returns
        -------
        dict
            A dictionary containing interpolated electron density (`n_e`),
            its time derivative (`dn_dt`), and the Greenwald fraction (`greenwald_fraction`).
        """
        if len(n_e) != len(t_n):
            raise CalculationError("n_e and t_n are different lengths")
        # get the gradient of n_E
        dn_dt = np.gradient(n_e, t_n)
        n_e = interp1(t_n, n_e, times)
        dn_dt = interp1(t_n, dn_dt, times)
        ip = -ip / 1e6  # Convert from A to MA and take positive value
        ip = interp1(t_ip, ip, times)
        a_minor = interp1(t_a, a_minor, times, bounds_error=False, fill_value=np.nan)
        # make sure aminor is not 0 or less than 0
        a_minor[a_minor <= 0] = 0.001
        n_g = abs(ip) / (np.pi * a_minor**2) * 1e20  # Greenwald density in m ^-3
        g_f = n_e / n_g
        return {"n_e": n_e, "dn_dt": dn_dt, "greenwald_fraction": g_f}

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
