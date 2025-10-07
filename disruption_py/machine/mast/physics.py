#!/usr/bin/env python3

"""
Physics methods for MAST.
"""

import numpy as np
from disruption_py.core.physics_method.errors import CalculationError
import xarray as xr
from disruption_py.core.utils.math import interp1
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.inout.xarray import XarrayConnection
from disruption_py.machine.tokamak import Tokamak


class MastPhysicsMethods:
    """
    MAST class to serve as a method holder.
    """

    @staticmethod
    @physics_method(
        columns=["ip"],
        tokamak=Tokamak.MAST,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """Get Ip parameters"""
        conn: XarrayConnection = params.mds_conn
        ip = conn.get_data(params.shot_id, "summary/ip")
        magtime = conn.get_data(params.shot_id, "summary/time")

        times = params.times
        ip = interp1(magtime, ip, times)
        return {"ip": ip}

    @staticmethod
    @physics_method(
        columns=["power_nbi", "power_radiated"],
        tokamak=Tokamak.MAST,
    )
    def get_power(params: PhysicsMethodParams):
        """Get power parameters"""
        conn: XarrayConnection = params.mds_conn
        power_nbi = conn.get_data(params.shot_id, "summary/power_nbi")
        power_radiated = conn.get_data(params.shot_id, "summary/power_radiated")
        base_time = conn.get_data(params.shot_id, "summary/time")

        times = params.times
        power_nbi = interp1(base_time, power_nbi, times)
        power_radiated = interp1(base_time, power_radiated, times)
        return {"power_nbi": power_nbi, "power_radiated": power_radiated}

    @staticmethod
    @physics_method(
        columns=["total_injected", "inboard_total", "outboard_total"],
        tokamak=Tokamak.MAST,
    )
    def get_gas(params: PhysicsMethodParams):
        """Get gas injection parameters"""
        conn: XarrayConnection = params.mds_conn
        total_injected = conn.get_data(params.shot_id, "gas_injection/total_injected")
        inboard_total = conn.get_data(params.shot_id, "gas_injection/inboard_total")
        outboard_total = conn.get_data(params.shot_id, "gas_injection/outboard_total")
        base_time = conn.get_data(params.shot_id, "gas_injection/time")

        times = params.times
        total_injected = interp1(base_time, total_injected, times)
        inboard_total = interp1(base_time, inboard_total, times)
        outboard_total = interp1(base_time, outboard_total, times)
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
        """Get Thomson parameters"""
        times = params.times
        conn: XarrayConnection = params.mds_conn
        t_e_core = conn.get_data(params.shot_id, "thomson_scattering/t_e_core")
        n_e_core = conn.get_data(params.shot_id, "thomson_scattering/n_e_core")
        base_time = conn.get_data(params.shot_id, "thomson_scattering/time")
        t_e_core = interp1(base_time, t_e_core, times)
        n_e_core = interp1(base_time, n_e_core, times)
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
            The parameters containing the MDSplus connection, shot id and more.

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
