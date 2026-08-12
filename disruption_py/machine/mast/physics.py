#!/usr/bin/env python3
# pylint: disable=duplicate-code

"""
Physics methods for MAST.
"""

import numpy as np
import scipy.constants as const

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.errors import (
    CalculationError,
    CustomError,
    DataError,
    MismatchCalculationError,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import causal_boxcar_smooth, gaussian_fit, interp1
from disruption_py.core.utils.misc import assert_equal_length
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
        magtime = params.get_data("summary/time")

        dip_dt = np.gradient(ip, magtime)

        times = params.times

        ip = MastUtilMethods.interpolate_1d(magtime, ip, times)
        dip_dt = MastUtilMethods.interpolate_1d(magtime, dip_dt, times)

        try:
            ip_prog = params.get_data("pulse_schedule/i_plasma", required=True)
            ip_prog_time = params.get_data("pulse_schedule/time", required=True)
        except DataError:
            # older shots lack a recorded pulse schedule; ip/dip_dt stand on their own
            ip_prog = [np.nan]
            dipprog_dt = [np.nan]
        else:
            dipprog_dt = np.gradient(ip_prog, ip_prog_time)
            ip_prog = MastUtilMethods.interpolate_1d(ip_prog_time, ip_prog, times)
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

        The Ohmic heating power is calculated from the dynamic loop voltage,
        inductive voltage, and plasma current.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing neutral beam injection power (`p_nbi`),
            total Ohmic heating power (`p_oh`), and radiated power (`p_rad`).
        """

        base_time = params.get_data("summary/time")
        times = params.times

        power_nbi = params.get_data("summary/power_nbi")
        power_nbi = MastUtilMethods.interpolate_1d(base_time, power_nbi, times)

        power_radiated = params.get_data("summary/power_radiated")
        power_radiated = MastUtilMethods.interpolate_1d(
            base_time, power_radiated, times
        )

        try:
            power_ohm = MastPhysicsMethods._get_p_ohm(params)
        except CustomError:
            # a missing equilibrium or Ip signal must not take out p_nbi and p_rad
            power_ohm = [np.nan]

        return {"p_nbi": power_nbi, "p_oh": power_ohm, "p_rad": power_radiated}

    @staticmethod
    def _get_p_ohm(params: PhysicsMethodParams):
        """
        Calculate the ohmic heating power from the dynamic loop voltage,
        inductive voltage, and plasma current.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        np.ndarray
            The total Ohmic heating power [W], on the requested time base.

        Raises
        ------
        CalculationError
            If any of the required signals are missing or misaligned.
        """
        # load relevant parameters
        r0 = params.get_data("equilibrium/magnetic_axis_r")
        li = params.get_data("equilibrium/li")
        v_loop = params.get_data("equilibrium/vloop_dynamic")
        ip = params.get_data("summary/ip")
        summary_time = params.get_data("summary/time")
        equilibrium_time = params.get_data("equilibrium/time")

        assert_equal_length(ip, summary_time)
        for signal in (r0, li, v_loop):
            assert_equal_length(signal, equilibrium_time)

        # compute derived quantities
        smooth_window_size = 30
        dip_dt = np.gradient(ip, summary_time)
        if len(dip_dt) >= smooth_window_size:
            dip_smoothed = causal_boxcar_smooth(dip_dt, smooth_window_size)
        else:
            dip_smoothed = dip_dt
        inductance = 4.0 * np.pi * 1.0e-7 * r0 * li / 2.0

        # interpolate to consistent time-base
        times = params.times
        v_loop = MastUtilMethods.interpolate_1d(equilibrium_time, v_loop, times)
        inductance = MastUtilMethods.interpolate_1d(equilibrium_time, inductance, times)
        dip_smoothed = MastUtilMethods.interpolate_1d(summary_time, dip_smoothed, times)
        ip = MastUtilMethods.interpolate_1d(summary_time, ip, times)

        # calculate p_oh
        v_inductive = inductance * dip_smoothed
        v_resistive = v_loop - v_inductive
        p_oh = ip * v_resistive

        # Set negative p_ohm values to 0
        return np.clip(p_oh, 0, None)

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

        times = params.times

        base_time = params.get_data("gas_injection/time")
        gas_total_injected = params.get_data("gas_injection/total_injected")
        gas_inboard_total = params.get_data("gas_injection/inboard_total")
        gas_outboard_total = params.get_data("gas_injection/outboard_total")

        gas_total_injected = MastUtilMethods.interpolate_1d(
            base_time, gas_total_injected, times
        )
        gas_inboard_total = MastUtilMethods.interpolate_1d(
            base_time, gas_inboard_total, times
        )
        gas_outboard_total = MastUtilMethods.interpolate_1d(
            base_time, gas_outboard_total, times
        )

        return {
            "gas_total_injected": gas_total_injected,
            "gas_inboard_total": gas_inboard_total,
            "gas_outboard_total": gas_outboard_total,
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

        base_time = params.get_data("thomson_scattering/time")
        te_core = params.get_data("thomson_scattering/t_e_core")
        ne_core = params.get_data("thomson_scattering/n_e_core")

        te_core = MastUtilMethods.interpolate_1d(base_time, te_core, times)
        ne_core = MastUtilMethods.interpolate_1d(base_time, ne_core, times)

        return {
            "te_core": te_core,
            "ne_core": ne_core,
        }

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

        n_e = params.get_data("summary/line_average_n_e", required=True)
        t_n = params.get_data("summary/time")
        ip = params.get_data("summary/ip")
        t_ip = params.get_data("summary/time")
        # a_minor/t_a: if equilibrium data is missing, the Greenwald fraction
        # below is returned as NaN while n_e/dn_dt still stand
        a_minor = params.get_data("equilibrium/minor_radius")
        t_a = params.get_data("equilibrium/time")
        times = params.times

        if len(n_e) != len(t_n):
            raise MismatchCalculationError(
                f"len(n_e) = {len(n_e)} vs. len(t_n) = {len(t_n)}"
            )
        assert_equal_length(n_e, t_n)
        assert_equal_length(ip, t_ip)
        # get the gradient of n_E
        dn_dt = np.gradient(n_e, t_n)
        n_e = interp1(t_n, n_e, times)
        dn_dt = interp1(t_n, dn_dt, times)
        ip = -ip / 1e6  # Convert from A to MA and take positive value
        ip = interp1(t_ip, ip, times)

        if not np.any(np.isfinite(a_minor)) or np.size(t_a) < 2:
            # No usable equilibrium: the densities stand, the Greenwald fraction cannot.
            params.logger.warning(
                "get_densities: equilibrium/minor_radius not available. "
                "Returning NaN for greenwald_fraction."
            )
            return {
                "n_e": n_e,
                "dn_dt": dn_dt,
                "greenwald_fraction": [np.nan],
            }

        assert_equal_length(a_minor, t_a)
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
        hcam = params.get_data(
            "soft_x_rays/horizontal_cam_upper",
            return_xarray=True,
            required=True,
        )

        hcam_channel = hcam.isel(horizontal_cam_upper_channel=0)
        hcam_channel = hcam_channel.squeeze(drop=True)
        sxr_time = hcam_channel.time.values
        sxr_core = hcam_channel.values

        times = params.times
        sxr_core = MastUtilMethods.interpolate_1d(sxr_time, sxr_core, times)

        hcam_channel = hcam.isel(horizontal_cam_upper_channel=7)
        hcam_channel = hcam_channel.squeeze(drop=True)
        sxr_edge = hcam_channel.values

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
            required=True,
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
        25% of the vertical scale from the magnetic axis). Only fan-array channels
        (channel_type == 0) are used.

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
        power = params.get_data("bolometer/power", required=True)
        rmag = params.get_data("equilibrium/magnetic_axis_r", required=True)

        bolo_time = params.get_data("bolometer/time")
        first_r = params.get_data("bolometer/first_point_r")
        first_z = params.get_data("bolometer/first_point_z")
        second_r = params.get_data("bolometer/second_point_r")
        second_z = params.get_data("bolometer/second_point_z")
        validity = params.get_data("bolometer/validity")
        channel_type = params.get_data("bolometer/channel_type")

        zmag = params.get_data("equilibrium/magnetic_axis_z")
        efit_time = params.get_data("equilibrium/time")

        assert_equal_length(power, bolo_time)
        assert_equal_length(rmag, efit_time)
        assert_equal_length(zmag, efit_time)

        if power.ndim != 2:
            raise MismatchCalculationError(
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
        # z_j(t) = fz + (Rmag(t) - fr) * (sz - fz) / (sr - fr)   [Rea et al. 2020, Eq. 2]
        d_r = sr - fr  # (n_fan,)
        d_z = sz - fz  # (n_fan,)

        with np.errstate(divide="ignore", invalid="ignore"):
            z_j = fz[:, np.newaxis] + (
                (rmag_t[np.newaxis, :] - fr[:, np.newaxis])
                * (d_z[:, np.newaxis] / d_r[:, np.newaxis])
            )  # (n_fan, n_times)

        dist_from_axis = np.abs(z_j - zmag_t[np.newaxis, :])  # (n_fan, n_times)

        valid_2d = fan_valid[:, np.newaxis] & ~np.isnan(z_j)
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

        prad_cva = MastUtilMethods.interpolate_1d(bolo_time, prad_cva, params.times)

        return {"prad_peaking": prad_cva}

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

        References
        -------
        - Rea et al. (2020), *Fusion Sci. Technol.* 76(8), 912-924.
          DOI: 10.1080/15361055.2020.1798589
        """
        # 2D profile arrays: dims (major_radius, time)
        te_xr = params.get_data(
            "thomson_scattering/t_e", return_xarray=True, required=True
        )
        ne_xr = params.get_data(
            "thomson_scattering/n_e", return_xarray=True, required=True
        )

        pe_xr = params.get_data(
            "thomson_scattering/p_e", return_xarray=True, required=True
        )

        te_profile = te_xr.values
        ne_profile = ne_xr.values
        if pe_xr is None:
            # p_e = n_e k T_e, with T_e in eV so that k T_e = e * T_e [J]
            pe_profile = ne_profile * te_profile * const.e
        else:
            pe_profile = pe_xr.values

        r_ts = te_xr.coords["major_radius"].values
        ts_time = te_xr.coords["time"].values

        rho, _, _ = MastUtilMethods.thomson_rho(params, r_ts)

        # channels within 30% of normalized radius (Rea et al. 2020, Eq. 7)
        core_mask = rho < 0.3
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

        te_pf = MastUtilMethods.interpolate_1d(ts_time, te_pf, params.times)
        ne_pf = MastUtilMethods.interpolate_1d(ts_time, ne_pf, params.times)
        pressure_pf = MastUtilMethods.interpolate_1d(ts_time, pressure_pf, params.times)

        return {
            "te_peaking": te_pf,
            "ne_peaking": ne_pf,
            "pressure_peaking": pressure_pf,
        }

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

        References
        -------
        - original source: [get_TS_data_cmod.m](https://github.com/MIT-PSFC/
        disruption-py/blob/matlab/CMOD/matlab-core/get_TS_data_cmod.m), adapted from
        the vertical C-Mod geometry to the radial MAST geometry.
        """

        # 2D profile array: dims (major_radius, time)
        te_xr = params.get_data("thomson_scattering/t_e", return_xarray=True)

        r_ts = te_xr.coords["major_radius"].values
        ts_time = te_xr.coords["time"].values
        te_profile = te_xr.values

        _, r_mag, a_minor = MastUtilMethods.thomson_rho(params, r_ts)

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
        te_hwhm = MastUtilMethods.interpolate_1d(ts_time, te_hwhm, params.times)
        return {"te_width": te_hwhm}

    @staticmethod
    @physics_method(
        columns=["z_error", "z_prog", "zcur", "v_z", "z_times_v_z"],
        tokamak=Tokamak.MAST,
    )
    def get_z_parameters(params: PhysicsMethodParams):
        """
        Retrieve the z parameters from the MAST diagnostic data.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters containing the Xarray connection, shot id and more.

        Returns
        -------
        dict
            A dictionary containing the z parameters: `z_error`,
            `z_prog`, `zcur`, `v_z`, and `z_times_v_z`.
        """
        z_ref = params.get_data("pulse_schedule/z_ref", required=True)
        zip_prx = params.get_data("controllers/zip_proxy", required=True)
        t_ctrl = params.get_data("controllers/time", required=True)
        ip_raw = params.get_data("summary/ip", required=True)
        t_ip = params.get_data("summary/time", required=True)

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
