#!/usr/bin/env python3
# pylint: disable=duplicate-code

"""
Physics methods for MAST.
"""

import numpy as np
import scipy.constants as const

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import (
    MismatchCalculationError,
    CalculationError,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import causal_boxcar_smooth, interp1
from disruption_py.core.utils.math import gaussian_fit, interp1
from disruption_py.inout.xr import XarrayDataConnection
from disruption_py.machine.mast.util import MastUtilMethods
from disruption_py.machine.tokamak import Tokamak


class MastPhysicsMethods:
    """
    This class provides methods to retrieve and calculate physics-related data
    for MAST.
    """

    @staticmethod
    @physics_method(
        columns=["ip", "dip_dt", "ip_prog", "dipprog_dt"],
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
        ip = params.get_data("summary/ip")
        ip_prog = params.get_data("pulse_schedule/i_plasma")
        ip_prog_time = params.get_data("pulse_schedule/time")
        magtime = params.get_data("summary/time")

        dip_dt = np.gradient(ip, magtime)
        dipprog_dt = np.gradient(ip_prog, ip_prog_time)

        times = params.times

        ip = MastUtilMethods.interpolate_1d(magtime, ip, times)
        ip_prog = MastUtilMethods.interpolate_1d(ip_prog_time, ip_prog, times)
        dip_dt = MastUtilMethods.interpolate_1d(magtime, dip_dt, times)
        dipprog_dt = MastUtilMethods.interpolate_1d(ip_prog_time, dipprog_dt, times)

        return {
            "ip": ip,
            "dip_dt": dip_dt,
            "ip_prog": ip_prog,
            "dipprog_dt": dipprog_dt,
        }

    @staticmethod
    @physics_method(
        columns=["p_nbi", "p_oh", "p_rad"],
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
            A dictionary containing neutral beam injection power (`p_nbi`),
            ohmic power (`p_oh`), and radiated power (`p_rad`).
        """

        power_nbi = params.get_data("summary/power_nbi")
        power_radiated = params.get_data("summary/power_radiated")
        base_time = params.get_data("summary/time")

        times = params.times
        power_nbi = MastUtilMethods.interpolate_1d(base_time, power_nbi, times)
        power_ohm = MastUtilMethods.interpolate_1d(base_time, power_ohm, times)
        power_radiated = MastUtilMethods.interpolate_1d(
            base_time, power_radiated, times
        )
        return {"p_nbi": power_nbi, "p_oh": power_ohm, "p_rad": power_radiated}

    @staticmethod
    @physics_method(
        columns=["gas_total_injected", "gas_inboard_total", "gas_outboard_total"],
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

        total_injected = params.get_data("gas_injection/total_injected")
        inboard_total = params.get_data("gas_injection/inboard_total")
        outboard_total = params.get_data("gas_injection/outboard_total")
        base_time = params.get_data("gas_injection/time")

        times = params.times
        total_injected = MastUtilMethods.interpolate_1d(
            base_time, total_injected, times
        )
        inboard_total = MastUtilMethods.interpolate_1d(base_time, inboard_total, times)
        outboard_total = MastUtilMethods.interpolate_1d(
            base_time, outboard_total, times
        )
        return {
            "gas_total_injected": total_injected,
            "gas_inboard_total": inboard_total,
            "gas_outboard_total": outboard_total,
        }

    @staticmethod
    @physics_method(
        columns=["te_core", "ne_core"],
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
            A dictionary containing core electron temperature (`te_core`) and
            core electron density (`ne_core`).
        """
        times = params.times

        t_e_core = params.get_data("thomson_scattering/t_e_core")
        n_e_core = params.get_data("thomson_scattering/n_e_core")
        base_time = params.get_data("thomson_scattering/time")

        te_core = MastUtilMethods.interpolate_1d(base_time, t_e_core, times)
        ne_core = MastUtilMethods.interpolate_1d(base_time, n_e_core, times)
        return {"te_core": te_core, "ne_core": ne_core}

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

        n_e = params.get_data("summary/line_average_n_e")
        t_n = params.get_data("summary/time")
        ip = params.get_data("summary/ip")
        t_ip = params.get_data("summary/time")
        a_minor = params.get_data("equilibrium/minor_radius")
        t_a = params.get_data("equilibrium/time")

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
            raise MismatchCalculationError(
                f"len(n_e) = {len(n_e)} vs. len(t_n) = {len(t_n)}"
            )
        MastUtilMethods.require_aligned("n_e", n_e, t_n)
        MastUtilMethods.require_aligned("ip", ip, t_ip)
        MastUtilMethods.require_aligned("a_minor", a_minor, t_a)
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
    @physics_method(
        columns=["sxr_core", "sxr_edge"],
        tokamak=Tokamak.MAST,
    )
    def get_sxr(params: PhysicsMethodParams):
        """
        Retrieve soft X-ray (SXR) data.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing SXR data (`sxr_data`) and
            corresponding time points (`sxr_time`).
        """
        hcam = params.get_data("soft_x_rays/horizontal_cam_upper", return_xarray=True)

        if hcam is not None:
            hcam_channel = hcam.isel(horizontal_cam_upper_channel=0)
            hcam_channel = hcam_channel.squeeze(drop=True)
            hcam_channel = hcam_channel.drop_vars(["horizontal_cam_upper_channel"])
            sxr_time = hcam_channel.time.values
            sxr_core = hcam_channel.values
        else:
            sxr_time = np.array([np.nan])
            sxr_core = np.array([np.nan])

        times = params.times
        sxr_core = MastUtilMethods.interpolate_1d(sxr_time, sxr_core, times)

        if hcam is not None:
            hcam_channel = hcam.isel(horizontal_cam_upper_channel=7)
            hcam_channel = hcam_channel.squeeze(drop=True)
            hcam_channel = hcam_channel.drop_vars(["horizontal_cam_upper_channel"])
            sxr_edge = hcam_channel.values
        else:
            sxr_edge = np.array([np.nan])

        times = params.times
        sxr_edge = MastUtilMethods.interpolate_1d(sxr_time, sxr_edge, times)

        return {"sxr_core": sxr_core, "sxr_edge": sxr_edge}

    @staticmethod
    @physics_method(
        columns=["d_alpha"],
        tokamak=Tokamak.MAST,
    )
    def get_dalpha(params: PhysicsMethodParams):
        """
        Retrieve D-alpha signal data.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing D-alpha signal data (`dalpha`).
        """

        dalpha = params.get_data(
            "spectrometer_visible/filter_spectrometer_dalpha_voltage",
            return_xarray=True,
        )

        dalpha = dalpha.isel(dalpha_channel=2)
        dalpha = dalpha.dropna(dim="time")
        dalpha = dalpha.squeeze(drop=True)
        dalpha = dalpha.drop_vars("dalpha_channel")

        dalpha_time = dalpha.time.values
        dalpha_data = dalpha.values

        times = params.times
        dalpha_data = MastUtilMethods.interpolate_1d(dalpha_time, dalpha_data, times)
        return {"d_alpha": dalpha_data}

    @staticmethod
    def _get_prad_peaking(
        times,
        power,
        bolo_time,
        first_r,
        first_z,
        second_r,
        second_z,
        validity,
        channel_type,
        rmag,
        zmag,
        efit_time,
    ):
        """
        Calculate the Prad CVA peaking factor from MAST bolometer data.

        Following Rea et al. (2020):
        - CVA: mean power in core channels / mean power in non-divertor channels

        Core channels are those whose chords pass within 6% of the vertical machine
        scale from the magnetic axis (Rea et al. 2020, Eq. 3). Divertor channels
        (excluded from the denominator) are those beyond 25% of the vertical scale.
        Only fan-array channels (channel_type == 0) are used.

        Parameters
        ----------
        times : array_like
            Requested time basis.
        power : np.ndarray
            Bolometer power array of shape (n_channels, n_times).
        bolo_time : np.ndarray
            Time base of the bolometer measurements.
        first_r, first_z : np.ndarray
            R and Z coordinates of the first point (detector/aperture) for each channel.
        second_r, second_z : np.ndarray
            R and Z coordinates of the second point (chord direction) for each channel.
        validity : np.ndarray
            Validity flags for each channel; 1 = valid, 0 = invalid.
        channel_type : np.ndarray
            Channel type identifier; 0 = vertical fan array.
        rmag, zmag : np.ndarray
            Magnetic axis R and Z positions on the equilibrium time base.
        efit_time : np.ndarray
            Time base of the equilibrium data.

        Returns
        -------
        dict
            A dictionary containing `prad_peaking`.
        """
        if power.ndim != 2:
            raise CalculationError(
                "Expected a 2-D (channel, time) bolometer power array, "
                f"got shape {power.shape}"
            )

        # Every per-channel array is indexed by the fan channels selected below
        n_channels = power.shape[0]
        per_channel = {
            "first_point_r": first_r,
            "first_point_z": first_z,
            "second_point_r": second_r,
            "second_point_z": second_z,
            "validity": validity,
            "channel_type": channel_type,
        }
        for name, array in per_channel.items():
            if array.shape != (n_channels,):
                raise CalculationError(
                    f"bolometer/{name} has shape {array.shape}, "
                    f"expected ({n_channels},) to match the power array"
                )

        # Use only fan-array channels (type 0) for consistent vertical coverage
        (fan_idx,) = np.where(channel_type == 0)
        if fan_idx.size == 0:
            raise CalculationError("No bolometer fan-array channels in this shot")

        fan_valid = validity[fan_idx] == 1  # (n_fan,)

        fr = first_r[fan_idx]
        fz = first_z[fan_idx]
        sr = second_r[fan_idx]
        sz = second_z[fan_idx]

        # Derive thresholds from fan geometry (Rea et al. 2020, Eq. 3)
        vert_scale = 2.0 * np.max(np.abs(sz))
        core_threshold = 0.06 * vert_scale
        div_threshold = 0.25 * vert_scale

        # Interpolate magnetic axis position onto the bolometer time base
        rmag_t = MastUtilMethods.interpolate_1d(efit_time, rmag, bolo_time)
        zmag_t = MastUtilMethods.interpolate_1d(efit_time, zmag, bolo_time)

        # Compute Z-intersection of each chord with vertical R = Rmag(t)
        # Z_j(t) = fz + (Rmag(t) - fr) * (sz - fz) / (sr - fr)   [Rea et al. 2020, Eq. 2]
        dR = sr - fr  # (n_fan,)
        dZ = sz - fz  # (n_fan,)

        with np.errstate(divide="ignore", invalid="ignore"):
            Z_j = fz[:, np.newaxis] + (
                (rmag_t[np.newaxis, :] - fr[:, np.newaxis])
                * (dZ[:, np.newaxis] / dR[:, np.newaxis])
            )  # (n_fan, n_times)

        dist_from_axis = np.abs(Z_j - zmag_t[np.newaxis, :])  # (n_fan, n_times)

        valid_2d = fan_valid[:, np.newaxis] & ~np.isnan(Z_j)
        core_mask = valid_2d & (dist_from_axis < core_threshold)
        all_but_div = valid_2d & ~(dist_from_axis > div_threshold)

        fan_power = power[fan_idx, :]  # (n_fan, n_times)

        core_vals = np.where(core_mask, fan_power, np.nan)
        abd_vals = np.where(all_but_div, fan_power, np.nan)
        core_count = np.isfinite(core_vals).sum(axis=0)
        abd_count = np.isfinite(abd_vals).sum(axis=0)
        core_mean = np.nansum(core_vals, axis=0) / np.where(
            core_count > 0, core_count, np.nan
        )
        abd_mean = np.nansum(abd_vals, axis=0) / np.where(
            abd_count > 0, abd_count, np.nan
        )

        prad_cva = np.where(abd_mean != 0, core_mean / abd_mean, np.nan)

        # Zero out times where equilibrium data is unavailable
        eq_nan = np.isnan(rmag_t) | np.isnan(zmag_t)
        prad_cva[eq_nan] = np.nan

        prad_cva = MastUtilMethods.interpolate_1d(bolo_time, prad_cva, times)

        return {"prad_peaking": prad_cva}

    @staticmethod
    @physics_method(
        columns=["prad_peaking"],
        tokamak=Tokamak.MAST,
    )
    def get_prad_peaking(params: PhysicsMethodParams):
        r"""
        Calculate the radiated power CVA peaking factor using the MAST bolometer arrays.

        Following Rea et al. (2020):

        $$
        P_{\text{rad,CVA}} = \frac{\langle P_j \rangle_{j \in C}}
                                   {\langle P_j \rangle_{j \notin D}}
        $$

        where C is the core bin (chords passing within 6% of the vertical machine scale
        from the magnetic axis) and D is the divertor bin (chords passing farther than
        25% of the vertical scale from the magnetic axis).

        The XDIV peaking factor is not computed for MAST because its symmetric
        double-null geometry makes a single divertor-vs-all metric ambiguous.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing `prad_peaking`.

        References
        -------
        - Rea et al. (2020), *Fusion Sci. Technol.* 76(8), 912-924.
          DOI: 10.1080/15361055.2020.1798589
        """
        conn: XarrayConnection = params.mds_conn
        shot_id = params.shot_id

        power = MastUtilMethods.require_signal(conn, shot_id, "bolometer/power")
        bolo_time = MastUtilMethods.require_time_base(conn, shot_id, "bolometer/time")
        first_r = MastUtilMethods.require_signal(
            conn, shot_id, "bolometer/first_point_r"
        )
        first_z = MastUtilMethods.require_signal(
            conn, shot_id, "bolometer/first_point_z"
        )
        second_r = MastUtilMethods.require_signal(
            conn, shot_id, "bolometer/second_point_r"
        )
        second_z = MastUtilMethods.require_signal(
            conn, shot_id, "bolometer/second_point_z"
        )
        validity = MastUtilMethods.require_signal(conn, shot_id, "bolometer/validity")
        channel_type = MastUtilMethods.require_signal(
            conn, shot_id, "bolometer/channel_type"
        )

        rmag = MastUtilMethods.require_signal(
            conn, shot_id, "equilibrium/magnetic_axis_r"
        )
        zmag = MastUtilMethods.require_signal(
            conn, shot_id, "equilibrium/magnetic_axis_z"
        )
        efit_time = MastUtilMethods.require_time_base(conn, shot_id, "equilibrium/time")

        MastUtilMethods.require_aligned("bolometer/power", power, bolo_time)
        MastUtilMethods.require_aligned("equilibrium/magnetic_axis_r", rmag, efit_time)
        MastUtilMethods.require_aligned("equilibrium/magnetic_axis_z", zmag, efit_time)

        return MastPhysicsMethods._get_prad_peaking(
            params.times,
            power,
            bolo_time,
            first_r,
            first_z,
            second_r,
            second_z,
            validity,
            channel_type,
            rmag,
            zmag,
            efit_time,
        )

    @staticmethod
    def _get_te_ne_peaking(times, te_profile, ne_profile, pe_profile, rho, ts_time):
        """
        Calculate Te, ne and pressure peaking factors from Thomson scattering
        profile data.

        Following Rea et al. (2020), the core bin is defined as channels with
        normalized effective radius rho < 0.3. The peaking factor is the ratio of
        the mean value in the core bin to the mean over all channels.

        Parameters
        ----------
        times : array_like
            Requested time basis.
        te_profile : np.ndarray
            Electron temperature profile, shape (n_channels, n_times), in eV.
        ne_profile : np.ndarray
            Electron density profile, shape (n_channels, n_times), in m^-3.
        pe_profile : np.ndarray
            Electron pressure profile, shape (n_channels, n_times), in Pa.
        rho : np.ndarray
            Normalized effective radius (rho = r/r_boundary) for each Thomson channel,
            shape (n_channels,).
        ts_time : np.ndarray
            Time base of the Thomson scattering measurements.

        Returns
        -------
        dict
            A dictionary containing `te_peaking`, `ne_peaking` and
            `pressure_peaking`.
        """
        core_mask = (
            rho < 0.3
        )  # channels within 30% of normalized radius (Rea et al. 2020, Eq. 7)
        n_times = len(ts_time)
        te_pf = np.full(n_times, np.nan)
        ne_pf = np.full(n_times, np.nan)
        pressure_pf = np.full(n_times, np.nan)

        for i_time in range(n_times):
            te_t = te_profile[:, i_time]
            ne_t = ne_profile[:, i_time]
            pe_t = pe_profile[:, i_time]

            # one shared mask, so the three factors are computed over the same channels
            valid = np.isfinite(te_t) & np.isfinite(ne_t) & (te_t > 0) & (ne_t > 0)
            core_valid = valid & core_mask

            if core_valid.sum() < 2 or valid.sum() < 3:
                continue

            te_avg = np.mean(te_t[valid])
            ne_avg = np.mean(ne_t[valid])
            pe_avg = np.mean(pe_t[valid])
            if te_avg > 0:
                te_pf[i_time] = np.mean(te_t[core_valid]) / te_avg
            if ne_avg > 0:
                ne_pf[i_time] = np.mean(ne_t[core_valid]) / ne_avg
            if pe_avg > 0:
                pressure_pf[i_time] = np.mean(pe_t[core_valid]) / pe_avg

        te_pf = MastUtilMethods.interpolate_1d(ts_time, te_pf, times)
        ne_pf = MastUtilMethods.interpolate_1d(ts_time, ne_pf, times)
        pressure_pf = MastUtilMethods.interpolate_1d(ts_time, pressure_pf, times)

        return {
            "te_peaking": te_pf,
            "ne_peaking": ne_pf,
            "pressure_peaking": pressure_pf,
        }

    @staticmethod
    @physics_method(
        columns=["te_peaking", "ne_peaking", "pressure_peaking"],
        tokamak=Tokamak.MAST,
    )
    def get_te_ne_peaking(params: PhysicsMethodParams):
        r"""
        Calculate Te, ne and pressure peaking factors from Thomson scattering
        profile data.

        The peaking factor for each quantity is defined following Rea et al. (2020):

        $$
        T_{e,\text{pf}} = \frac{\langle T_j \rangle_{j \in C}}{T_{\text{avg}}}
        \qquad
        n_{e,\text{pf}} = \frac{\langle n_j \rangle_{j \in C}}{n_{\text{avg}}}
        \qquad
        p_{e,\text{pf}} = \frac{\langle p_j \rangle_{j \in C}}{p_{\text{avg}}}
        $$

        where C is the set of Thomson scattering channels with normalized effective
        radius $\rho_j < 0.3$ and the denominators are mean values over all channels.

        The electron pressure profile is read from ``thomson_scattering/p_e`` where
        available, and otherwise reconstructed as $p_e = n_e k T_e$.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing `te_peaking`, `ne_peaking` and
            `pressure_peaking`.

        Raises
        ------
        CalculationError
            If Thomson scattering profile data is not available for this shot.
            Required zarr paths: ``thomson_scattering/t_e``,
            ``thomson_scattering/n_e``, ``equilibrium/magnetic_axis_r``,
            ``equilibrium/minor_radius``.

        References
        -------
        - Rea et al. (2020), *Fusion Sci. Technol.* 76(8), 912-924.
          DOI: 10.1080/15361055.2020.1798589
        """
        conn: XarrayConnection = params.mds_conn

        # 2D profile arrays: dims (major_radius, time)
        te_xr = conn.get_data(
            params.shot_id, "thomson_scattering/t_e", return_xarray=True
        )
        ne_xr = conn.get_data(
            params.shot_id, "thomson_scattering/n_e", return_xarray=True
        )
        pe_xr = conn.get_data(
            params.shot_id, "thomson_scattering/p_e", return_xarray=True
        )

        if any(x is None for x in (te_xr, ne_xr)):
            raise CalculationError(
                "Thomson scattering profile data not available. "
                "Requires zarr paths: thomson_scattering/t_e, thomson_scattering/n_e."
            )

        te_profile = te_xr.values.squeeze()
        ne_profile = ne_xr.values.squeeze()
        if pe_xr is None:
            # p_e = n_e k T_e, with T_e in eV so that k T_e = e * T_e [J]
            params.logger.warning(
                "thomson_scattering/p_e unavailable, reconstructing from n_e and t_e."
            )
            pe_profile = ne_profile * te_profile * const.e
        else:
            pe_profile = pe_xr.values.squeeze()

        r_ts = te_xr.coords["major_radius"].values
        ts_time = te_xr.coords["time"].values

        rho, _, _ = MastUtilMethods.thomson_rho(conn, params.shot_id, r_ts)

        return MastPhysicsMethods._get_te_ne_peaking(
            params.times,
            te_profile,
            ne_profile,
            pe_profile,
            rho,
            ts_time,
        )

    @staticmethod
    def _get_te_width(times, te_profile, r_ts, ts_time, r_mag, a_minor):
        """
        Fit a Gaussian to each electron temperature profile and return its
        half-width at half-maximum.

        Parameters
        ----------
        times : array_like
            Requested time basis.
        te_profile : np.ndarray
            Electron temperature profile, shape (n_channels, n_times), in eV.
        r_ts : np.ndarray
            Major radius of each Thomson scattering channel [m], shape (n_channels,).
        ts_time : np.ndarray
            Time base of the Thomson scattering measurements.
        r_mag : float
            Time-averaged major radius of the magnetic axis [m], used to reject fits
            whose centre falls outside the plasma.
        a_minor : float
            Time-averaged minor radius [m], used as the scale for both the centre and
            the width rejection windows.

        Returns
        -------
        dict
            A dictionary containing the electron temperature profile width
            (`te_width`).
        """
        # sort by major radius
        idx = np.argsort(r_ts)
        r_ts = r_ts[idx]
        te_profile = te_profile[idx]
        # init output
        te_hwhm = np.full(len(ts_time), np.nan)
        # select valid times, dropping the pre-shot part of the TS time base
        (valid_times,) = np.where(ts_time > 0)
        for i_time in valid_times:
            y = te_profile[:, i_time]
            (ok_indices,) = np.where(np.isfinite(y) & (y > 0))
            # skip if not enough points
            if len(ok_indices) < 3:
                continue
            # working arrays
            y = y[ok_indices]
            r = r_ts[ok_indices]
            # initial guess
            i = y.argmax()
            guess = [y[i], r[i], (r.max() - r.min()) / 3]
            # actual fit; MAST level-2 data carries no Te uncertainty, so unlike the
            # C-Mod equivalent the fit is unweighted
            try:
                _, pmean, psigma = gaussian_fit(r, y, guess)
            except RuntimeError as exc:
                if str(exc).startswith("Optimal parameters not found"):
                    continue
                raise exc
            # reject points whose fitted centre falls outside the plasma
            if np.abs(pmean - r_mag) > a_minor:
                continue
            # store output
            te_hwhm[i_time] = np.abs(psigma)
        # rescale from sigma to HWHM
        # https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        te_hwhm *= np.sqrt(2 * np.log(2))
        # reject points with unphysical HWHM
        te_hwhm[te_hwhm > a_minor] = np.nan
        # time interpolation
        te_hwhm = MastUtilMethods.interpolate_1d(ts_time, te_hwhm, times)
        return {"te_width": te_hwhm}

    @staticmethod
    @physics_method(columns=["te_width"], tokamak=Tokamak.MAST)
    def get_te_width(params: PhysicsMethodParams):
        """
        Retrieve the electron temperature profile from the Thomson scattering (TS)
        diagnostic, then calculate the half-width at half-maximum of the Gaussian
        fit to the profile.

        MAST's TS channels lie along the major radius, so the fit is performed
        against `major_radius` rather than the vertical coordinate used on C-Mod.
        Fits are discarded when the fitted centre falls further than a minor radius
        from the magnetic axis, or when the resulting width exceeds a minor radius.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the electron temperature profile width
            (`te_width`).

        Raises
        ------
        CalculationError
            If Thomson scattering profile data is not available for this shot.
            Required zarr paths: ``thomson_scattering/t_e``,
            ``equilibrium/magnetic_axis_r``, ``equilibrium/minor_radius``.

        References
        -------
        - original source: [get_TS_data_cmod.m](https://github.com/MIT-PSFC/
        disruption-py/blob/matlab/CMOD/matlab-core/get_TS_data_cmod.m), adapted from
        the vertical C-Mod geometry to the radial MAST geometry.
        """
        conn: XarrayConnection = params.mds_conn

        # 2D profile array: dims (major_radius, time)
        te_xr = conn.get_data(
            params.shot_id, "thomson_scattering/t_e", return_xarray=True
        )
        if te_xr is None:
            raise CalculationError(
                "Thomson scattering profile data not available. "
                "Requires zarr path: thomson_scattering/t_e."
            )

        r_ts = te_xr.coords["major_radius"].values
        ts_time = te_xr.coords["time"].values

        _, r_mag, a_minor = MastUtilMethods.thomson_rho(conn, params.shot_id, r_ts)

        return MastPhysicsMethods._get_te_width(
            params.times,
            te_xr.values.squeeze(),
            r_ts,
            ts_time,
            r_mag,
            a_minor,
        )

    @staticmethod
    @physics_method(
        columns=["z_error", "z_prog", "zcur", "v_z", "z_times_v_z"],
        tokamak=Tokamak.MAST,
    )
    def get_z_parameters(params: PhysicsMethodParams):
        conn: XarrayConnection = params.mds_conn

        z_ref = conn.get_data(params.shot_id, "pulse_schedule/z_ref").squeeze()
        zip_prx = conn.get_data(params.shot_id, "controllers/zip_proxy").squeeze()
        t_ctrl = conn.get_data(params.shot_id, "controllers/time")
        ip_raw = conn.get_data(params.shot_id, "summary/ip").squeeze()
        t_ip = conn.get_data(params.shot_id, "summary/time")

        if any(
            not np.isfinite(x).any() for x in (z_ref, zip_prx, t_ctrl, ip_raw, t_ip)
        ):
            params.logger.warning(
                "z_parameters: required signal(s) not available. Returning NaNs."
            )
            return {
                col: [np.nan]
                for col in ("z_error", "z_prog", "zcur", "v_z", "z_times_v_z")
            }

        ip_ctrl = MastUtilMethods.interpolate_1d(t_ip, ip_raw, t_ctrl)

        # Avoid amplifying noise when plasma is off (threshold: 10 kA)
        safe_ip = np.where(np.abs(ip_ctrl) > 1e4, ip_ctrl, np.nan)

        zcur = zip_prx / safe_ip  # [m·A / A = m]
        z_error = zcur - z_ref
        v_z = np.gradient(zcur, t_ctrl)
        z_times_v_z = zcur * v_z

        return {
            "z_prog": MastUtilMethods.interpolate_1d(t_ctrl, z_ref, params.times),
            "zcur": MastUtilMethods.interpolate_1d(t_ctrl, zcur, params.times),
            "z_error": MastUtilMethods.interpolate_1d(t_ctrl, z_error, params.times),
            "v_z": MastUtilMethods.interpolate_1d(t_ctrl, v_z, params.times),
            "z_times_v_z": MastUtilMethods.interpolate_1d(
                t_ctrl, z_times_v_z, params.times
            ),
        }
