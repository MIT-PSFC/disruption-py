from disruption_py.shots.shot import Shot, DEFAULT_SHOT_COLUMNS
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources
import logging
import traceback

import pandas as pd
import numpy as np
import scipy
import netCDF4 as nc

import MDSplus
from MDSplus import *

from disruption_py.utils import interp1, interp2, smooth, gaussian_fit, gsastd, get_bolo, power, efit_rz_interp
import disruption_py.data
D3D_DISRUPTED_SHOT = 175552
# Retrieve efit from EFIT01, Peaking Factor Nodes: dpsrdcva dpsrdxdiv dpstepf dpsnepf
D3D_PEAKING_FACTORS_SHOT = 180808
D3D_MAX_SHOT_TIME = 7.0  # [s]
DEFAULT_A_MINOR = 0.56  # [m]
"""
Useful Examples:
https://diii-d.gat.com/diii-d/MDSplusAPI_Python_pg1
https://diii-d.gat.com/diii-d/Gadata_py
"""


class D3DShot(Shot):
    # Tokamak Variables
    nominal_flattop_radius = 0.59
    # EFIT Variables
    efit_cols = {'beta_n': r'\efit_a_eqdsk:betan', 'beta_p': r'\efit_a_eqdsk:betap', 'kappa': r'\efit_a_eqdsk:kappa', 'li': r'\efit_a_eqdsk:li', 'upper_gap': r'\efit_a_eqdsk:gaptop', 'lower_gap': r'\efit_a_eqdsk:gapbot',
                 'q0': r'\efit_a_eqdsk:q0', 'qstar': r'\efit_a_eqdsk:qstar', 'q95': r'\efit_a_eqdsk:q95', 'Wmhd': r'\efit_a_eqdsk:wmhd', 'chisq': r'\efit_a_eqdsk:chisq'}
    # 'v_loop_efit': ,r'\efit_a_eqdsk:vsurf', 'bt0': r'\efit_a_eqdsk:bt0'
    efit_derivs = {'beta_p': 'dbetap_dt', 'li': 'dli_dt', 'Wmhd': 'dWmhd_dt'}
    rt_efit_cols = {'beta_p_RT': r'\efit_a_eqdsk:betap', 'li_RT': r'\efit_a_eqdsk:li',
                    'q95_RT': r'\efit_a_eqdsk:q95', 'Wmhd_RT': r'\efit_a_eqdsk:wmhd', 'chisq': r'\efit_a_eqdsk:chisq'}
    # 'v_loop_efit_RT': r'\efit_a_eqdsk:vsurf',

    # Disruption Variables
    dt_before_disruption = 0.002
    duration_before_disruption = 0.10

    def __init__(self, shot_id, efit_tree_name, data=None, times=None, disruption_time=None, override_cols=True, **kwargs):
        super().__init__(shot_id, data)
        self._times = times
        self.conn = MDSplus.Connection('atlas.gat.com')
        self.efit_tree_name = str(efit_tree_name)
        self.disruption_time = disruption_time
        self.disrupted = self.disruption_time is not None
        self.override_cols = override_cols
        self.data = data
        timebase_signal = kwargs.pop('timebase_signal', None)
        populate = kwargs.pop('populate', 'default')
        if self.data is not None and self._times is None:
            try:
                self._times = self.data['time'].to_numpy()
                # Check if the timebase is in ms instead of s
                if self._times[-1] > D3D_MAX_SHOT_TIME:
                    self._times /= 1000  # [ms] -> [s]
            except KeyError as e:
                self.logger.warning(
                    f"[Shot {self._shot_id}]:Shot constructor was passed data but no timebase.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]:{traceback.format_exc()}")
                self._times = self.get_timebase(timebase_signal, **kwargs)
        if self._times is None:
            self._times = self.get_timebase(timebase_signal, **kwargs)
        if populate == 'l_mode':
            self._populate_l_mode_data()
        else:
            self._populate_shot_data(data is not None)

    def _populate_l_mode_data(self):
        self.data = pd.concat([self.get_efit_parameters(), self.get_density_parameters(
        ), self.get_peaking_factors(), self.get_kappa_area(), self.get_time_until_disrupt()], axis=1)
        self.data['time'] = self._times
        self.data['shot'] = self._shot_id

    def _populate_shot_data(self, already_populated=False):
        local_data = pd.concat([self.get_efit_parameters(), self.get_rt_efit_parameters(), self.get_density_parameters(), self.get_rt_density_parameters(), self.get_ip_parameters(
        ), self.get_rt_ip_parameters(), self.get_power_parameters(), self.get_z_parameters(), self.get_zeff_parameters(), self.get_shape_parameters(), self.get_kappa_area(), self.get_time_until_disrupt()], axis=1)
        local_data = local_data.loc[:, ~local_data.columns.duplicated()]
        if not already_populated:
            self.data = local_data
            self.data['time'] = self._times
            self.data['shot'] = self._shot_id

    def get_timebase(self, timebase_signal, **kwargs):
        if timebase_signal == None or timebase_signal == 'disruption_timebase':
            minimum_ip = kwargs.get('minimum_ip', 400.e3)
            minimum_duration = kwargs.get('minimum_duration', 0.1)
            return self.get_disruption_timebase(minimum_ip, minimum_duration)
        elif timebase_signal == 'ip':
            _, ip_time = self._get_signal(
                f"ptdata('ip', {self._shot_id})", interpolate=False)
            return ip_time
        elif timebase_signal == 'flattop':
            _, ip_time = self._get_signal(
                f"ptdata('ip', {self._shot_id})", interpolate=False)
            return self.get_flattop_timebase(ip_time)
        else:
            raise NotImplementedError(
                "Only 'disruption_timebase' and 'ip' are supported for timebase_signal.")

    def get_disruption_timebase(self, minimum_ip=400.e3, minimum_duration=0.1):
        self.conn.openTree('d3d', self._shot_id)
        raw_ip, ip_time = self._get_signal(
            f"ptdata('ip', {self._shot_id})", interpolate=False)
        baseline = np.mean(raw_ip[0:10])
        ip = raw_ip - baseline
        duration, ip_max = self.get_end_of_shot(ip, ip_time, 100e3)
        if duration < minimum_duration or np.abs(ip_max) < minimum_ip:
            raise NotImplementedError()
        times = np.arange(0.100, duration+0.025, 0.025)
        if self.disrupted:
            additional_times = np.arange(
                self.disruption_time-self.duration_before_disruption, self.disruption_time + self.dt_before_disruption, self.dt_before_disruption)
            times = times[np.where(times < (self.disruption_time -
                                            self.duration_before_disruption))]
            times = np.concatenate((times, additional_times))
        else:
            ip_start = np.argmax(ip_time <= .1)
            ip_end = np.argmax(raw_ip[ip_start:] <= 100000) + ip_start
            return ip_time[ip_start:ip_end]  # [ms] -> [s]
        return times

    def get_end_of_shot(self, signal, signal_time, threshold=1.e5):
        duration = 0
        signal_max = 0
        if threshold < 0:
            raise Warning("Threshold is negative.")
        base_indices = np.where(signal_time <= 0.0)
        if len(base_indices) > 0:
            baseline = np.mean(signal[base_indices])
        else:
            baseline = 0
        signal = signal - baseline
        # Check if there was a finite signal otherwise consider the shot a "no plasma" shot
        finite_indices = np.where(
            (signal_time >= 0.0) & (np.abs(signal) > threshold))
        if len(finite_indices) == 0:
            return duration, signal_max
        else:
            dt = np.diff(signal_time)
            duration = np.sum(dt[finite_indices[:-1]])
            if duration < 0.1:  # Assuming < 100 ms is not a bona fide plasma
                duration = 0
                return duration, signal_max
        polarity = np.sign(
            np.trapz(signal[finite_indices], signal_time[finite_indices]))
        polarized_signal = polarity * signal
        valid_indices = np.where(
            (polarized_signal >= threshold) & (signal_time > 0.0))
        duration = signal_time[np.max(valid_indices)]
        if len(valid_indices) == signal_time.size:
            duration = - duration
        signal_max = np.max(polarized_signal)*polarity
        return duration, signal_max

    def get_flattop_timebase(self, times):
        # time_until_disrupt = self.disruption_time - times
        # time_until_disrupt[time_until_disrupt < 0] = np.nan
        # indices_no_disrupt = np.where(np.isnan(time_until_disrupt))
        # indices_disrupt = np.where(~np.isnan(time_until_disrupt))
        try:
            t_ip_prog = self.conn.get(
                f"dim_of(ptdata('iptipp', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_prog = self.conn.get(
                f"ptdata('iptipp', {self._shot_id})").data()  # [A]
            polarity = np.unique(self.conn.get(
                f"ptdata('iptdirect', {self._shot_id})").data())
            if len(polarity) > 1:
                self.logger.info(
                    f"[Shot {self._shot_id}]:Polarity of Ip target is not constant. Using value at first timestep.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]: Polarity array {polarity}")
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, times, 'linear')
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get programmed plasma current parameters")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        epsoff = self.conn.get(f"ptdata('epsoff', {self._shot_id})").data()
        t_epsoff = self.conn.get(
            f"dim_of(ptdata('epsoff', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
        t_epsoff += .001  # Avoid problem with simultaneity of epsoff being triggered exactly on the last time sample
        epsoff = interp1(t_epsoff, epsoff, times, 'linear')
        railed_indices = np.where(np.abs(epsoff) > .5)
        power_supply_railed = np.zeros(len(times))
        power_supply_railed[railed_indices] = 1
        indices_flattop = np.where((np.abs(dipprog_dt) <= 2.e3) & (
            np.abs(ip_prog) > 100e3) & (power_supply_railed != 1))
        return times[indices_flattop]

    def get_time_until_disrupt(self):
        if self.disrupted:
            return pd.DataFrame({'time_until_disrupt': self.disruption_time - self._times})
        return pd.DataFrame({'time_until_disrupt': np.full(self._times.size, np.nan)})

    def get_efit_parameters(self):
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        efit_data = {k: self.conn.get(v).data()
                     for k, v in self.efit_cols.items()}
        efit_time = self.conn.get(
            r'\efit_a_eqdsk:atime').data()/1.e3  # [ms] -> [s]
        if self._times is None:
            self._times = efit_time
        # EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of EFIT parameters that can indicate
        # invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        invalid_indices = np.where(efit_data['chisq'] > 50)
        del efit_data['chisq']
        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        for param in self.efit_derivs:
            efit_data[self.efit_derivs[param]] = np.gradient(
                efit_data[param], efit_time)
        if not np.array_equal(self._times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(
                    efit_time, efit_data[param], self._times)
        return pd.DataFrame(efit_data)

    def get_rt_efit_parameters(self):
        self.conn.openTree('efitrt1', self._shot_id)
        efit_data = {k: self.conn.get(v).data()
                     for k, v in self.rt_efit_cols.items()}
        efit_time = self.conn.get(
            r'\efit_a_eqdsk:atime').data()/1.e3  # [ms] -> [s]
        # EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of EFIT parameters that can indicate
        # invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        invalid_indices = np.where(efit_data['chisq'] > 50)
        del efit_data['chisq']
        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        if not np.array_equal(self._times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(
                    efit_time, efit_data[param], self._times)
        return pd.DataFrame(efit_data)

    def get_H_parameters(self):
        self.conn.openTree('transport', self._shot_id)
        try:
            h_98, _ = self._get_signal(r'\H_THH98Y2')
        except ValueError as e:
            self.logger.info(
                f"[Shot {self._shot_id}]: Failed to get H98 signal. Returning NaNs.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            h_98 = np.full(self._times.size, np.nan)
        self.conn.openTree('d3d', self._shot_id)
        try:
            h_alpha, _ = self._get_signal(r'\fs04')
        except ValueError as e:
            self.logger.info(
                f"[Shot {self._shot_id}]: Failed to get H_alpha signal. Returning NaNs.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            h_alpha = np.full(self._times.size, np.nan)
        return pd.DataFrame({'H98': h_98, 'H_alpha': h_alpha})

    def get_power_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        # Get neutral beam injected power
        try:
            p_nbi, t_nbi = self._get_signal(
                r'\d3d::top.nb:pinj', interpolate=False)
            p_nbi = p_nbi.astype(np.float64)
            p_nbi *= 1.e3  # [KW] -> [W]
            if len(t_nbi) > 2:
                p_nbi = interp1(t_nbi, p_nbi, self._times,
                                'linear', bounds_error=False, fill_value=0.)
            else:
                self.logger.info(
                    f"[Shot {self._shot_id}]:No NBI power data found in this shot.")
                p_nbi = np.zeros(len(self._times))
        except MdsException as e:
            p_nbi = np.zeros(len(self._times))
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to open NBI node")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Get electron cycholotrn heating (ECH) power. It's poitn data, so it's not stored in an MDSplus tree
        self.conn.openTree('rf', self._shot_id)
        try:
            p_ech, t_ech = self._get_signal(
                r'\top.ech.total:echpwrc', interpolate=False)
            if len(t_ech) > 2:
                p_ech = interp1(t_ech, p_ech, self._times,
                                'linear', bounds_error=False, fill_value=0.)
            else:
                self.logger.info(
                    f"[Shot {self._shot_id}]:No ECH power data found in this shot. Setting to zeros")
                p_ech = np.zeros(len(self._times))
        except MdsException as e:
            p_ech = np.zeros(len(self._times))
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to open ECH node. Setting to zeros")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Get ohmic power and loop voltage
        p_ohm, v_loop = self.get_ohmic_parameters()
        # Radiated power
        # We had planned to use the standard signal r'\bolom::prad_tot' for this
        # parameter.  However, the processing involved in calculating \prad_tot
        # from the arrays of bolometry channels involves non-causal filtering with
        # a 50 ms window.  This is not acceptable for our purposes.  Tony Leonard
        # provided us with the two IDL routines that are used to do the automatic
        # processing that generates the \prad_tot signal in the tree (getbolo.pro
        # and powers.pro).  I converted them into Matlab routines, and modified the
        # analysis so that the smoothing is causal, and uses a shorter window.
        smoothing_window = 0.010  # [s]
        try:
            self.conn.openTree("bolom", self._shot_id)
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to open bolom tree.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        bol_prm, _ = self._get_signal(r'\bol_prm', interpolate=False)
        lower_channels = [f"bol_u{i+1:02d}_v" for i in range(24)]
        upper_channels = [f"bol_l{i+1:02d}_v" for i in range(24)]
        bol_channels = lower_channels + upper_channels
        bol_signals = []
        bol_times = []
        for i in range(48):
            bol_signal, bol_time = self._get_signal(
                fr"\top.raw:{bol_channels[i]}", interpolate=False)
            bol_signals.append(bol_signal)
            bol_times.append(bol_time)
        a_struct = get_bolo(self._shot_id, bol_channels,
                            bol_prm, bol_signals, bol_times)
        ier = 0
        for j in range(48):
            # TODO: Ask about how many valid channels are needed for proper calculation
            if a_struct.channels[j].ier == 1:
                ier = 1
                p_rad = np.full(len(self._times), np.nan)
                break
        if ier == 0:
            b_struct = power(a_struct)
            p_rad = b_struct.pwrmix  # [W]
            p_rad = interp1(a_struct.time, p_rad, self._times, 'linear')

        # Remove any negative values from the power data
        p_rad[np.isinf(p_rad)] = np.nan
        p_rad[p_rad < 0] = 0
        p_nbi[p_nbi < 0] = 0
        p_ech[p_ech < 0] = 0

        p_input = p_rad + p_nbi + p_ech  # [W]
        rad_fraction = p_rad/p_input
        rad_fraction[np.isinf(rad_fraction)] = np.nan

        # Computer P_sol, defined as P_in - P_rad
        p_sol = p_input - p_rad

        return pd.DataFrame({'p_rad': p_rad, 'p_nbi': p_nbi, 'p_ech': p_ech, 'p_ohm': p_ohm, 'radiated_fraction': rad_fraction, 'v_loop': v_loop})

    def get_ohmic_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        # Get edge loop voltage and smooth it a bit with a median filter
        try:
            v_loop, t_v_loop = self._get_signal(
                f'ptdata("vloopb", {self._shot_id})', interpolate=False)
            v_loop = scipy.signal.medfilt(v_loop, 11)
            v_loop = interp1(t_v_loop, v_loop, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to open VLOOPB node. Setting to NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            v_loop = np.full(len(self._times), np.nan)
            t_v_loop = v_loop.copy()
       # Get plasma current
        try:
            ip = self.conn.get(f"ptdata('ip', {self._shot_id})").data()  # [A]
            t_ip = self.conn.get(
                f"dim_of(ptdata('ip', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            # We choose a 20-point width for gsastd. This means a 10ms window for ip smoothing
            dipdt_smoothed = gsastd(t_ip, ip, 1, 20, 3, 1, 0)
            self.conn.openTree(self.efit_tree_name, self._shot_id)
            li, t_li = self._get_signal(r'\efit_a_eqdsk:li', interpolate=False)
            chisq = self.conn.get(r'\efit_a_eqdsk:chisq').data()
            # Filter out invalid indices of efit reconstruction
            invalid_indices = None  # TODO: Finish
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Unable to get plasma current data. p_ohm set to NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            p_ohm = np.full(len(self._times), np.nan)
            return pd.DataFrame({'p_ohm': p_ohm, 'v_loop': v_loop})
        # [m] For simplicity, use fixed r_0 = 1.67 for DIII-D major radius
        r_0 = 1.67
        inductance = 4.*np.pi*r_0 * li/2  # [H]
        inductance = interp1(t_li, inductance, self._times, 'linear')
        ip = interp1(t_ip, ip, self._times, 'linear')
        dipdt_smoothed = interp1(t_ip, dipdt_smoothed, self._times, 'linear')
        v_inductive = inductance * dipdt_smoothed  # [V]
        v_resistive = v_loop - v_inductive  # [V]
        p_ohm = ip * v_resistive  # [W]
        return pd.DataFrame({'p_ohm': p_ohm, 'v_loop': v_loop})

    def get_density_parameters(self):
        ne = np.full(len(self._times), np.nan)
        g_f = ne.copy()
        dne_dt = ne.copy()
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        try:
            t_ne = self.conn.get(
                r'dim_of(\density)').data()/1.e3  # [ms] -> [s]
            ne = self.conn.get(r'\density').data()*1.e6  # [cm^3] -> [m^3]
            dne_dt = np.gradient(ne, t_ne)
            # NOTE: t_ne has higher resolution than efit_time so t_ne[0] < efit_time[0] because of rounding, meaning we need to allow extrapolation
            ne = interp1(t_ne, ne, self._times, 'linear',
                         bounds_error=False, fill_value='extrapolate')
            dne_dt = interp1(t_ne, dne_dt, self._times, 'linear',
                             bounds_error=False, fill_value='extrapolate')
            t_ip = self.conn.get(
                f"dim_of(ptdata('ip', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip = self.conn.get(f"ptdata('ip', {self._shot_id})").data()  # [A]
            ipsign = np.sign(np.sum(ip))
            ip = interp1(t_ip, ip*ipsign, self._times, 'linear')
            a_minor = self.conn.get(r'\efit_a_eqdsk:aminor').data()  # [m]
            t_a = self.conn.get(
                r'\efit_a_eqdsk:atime').data()/1.e3  # [ms] -> [s]
            a_minor = interp1(t_a, a_minor, self._times, 'linear')
            with np.errstate(divide='ignore'):
                n_g = ip/1.e6 / (np.pi*a_minor**2)  # [MA/m^2]
                g_f = ne/1.e20 / n_g  # TODO: Fill in units
        except MdsException as e:
            # TODO: Confirm that there is a separate exception if ptdata name doesn't exist
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get some parameter")
            self.logger.debug(
                f"[Shot {self._shot_id}]::{traceback.format_exc()}")
        return pd.DataFrame({'n_e': ne, 'Greenwald_fraction': g_f, 'dn_dt': dne_dt})

    def get_rt_density_parameters(self):
        ne_rt = np.full(len(self._times), np.nan)
        g_f_rt = ne_rt.copy()
        dne_dt_rt = ne_rt.copy()
        try:
            t_ne_rt = self.conn.get(
                f"dim_of(ptdata('dssdenest', {self._shot_id}))").data()/1.e3  # [ms] to [s]
            # [10^19 m^-3] -> [m^-3]
            ne_rt = self.conn.get(
                f"ptdata('dssdenest', {self._shot_id})").data()*1.e19
            dne_dt_rt = np.gradient(ne_rt, t_ne_rt)  # [m^-3/s]
            ne_rt = interp1(t_ne_rt, ne_rt, self._times, 'linear')
            dne_dt_rt = interp1(t_ne_rt, dne_dt_rt, self._times, 'linear')
            try:
                t_ip_rt = self.conn.get(
                    f"dim_of(ptdata('ipsip', {self._shot_id}))").data()/1.e3  # [ms] to [s]
                ip_rt = self.conn.get(
                    f"ptdata('ipsip', {self._shot_id})").data()  # [MA]
            except Exception as e:
                t_ip_rt = self.conn.get(
                    f"dim_of(ptdata('ipspr15v', {self._shot_id}))").data()/1.e3  # [ms] to [s]
                ip_rt = self.conn.get(
                    f"ptdata('ipspr15v', {self._shot_id})").data()  # [MA]
            ip_sign = np.sign(np.sum(ip_rt))
            ip = interp1(t_ip_rt, ip_rt*ip_sign, self._times, 'linear')
            self.conn.openTree('efitrt1', self._shot_id)
            a_minor_rt = self.conn.get(r'\efit_a_eqdsk:aminor').data()  # [m]
            t_a_rt = self.conn.get(
                r'\efit_a_eqdsk:atime').data()/1.e3  # [ms] -> [s]
            a_minor_rt = interp1(t_a_rt, a_minor_rt, self._times, 'linear')
            with np.errstate(divide='ignore'):
                n_g_rt = ip/1.e6 / (np.pi*a_minor_rt**2)  # [MA/m^2]
                g_f_rt = ne_rt/1.e20 / n_g_rt  # TODO: Fill in units
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get some parameter")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # ' dne_dt_RT': dne_dt_rt
        return pd.DataFrame({'n_e_RT': ne_rt, 'Greenwald_fraction_RT': g_f_rt})

    def get_ip_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        ip = np.full(len(self._times), np.nan)
        ip_prog = np.full(len(self._times), np.nan)
        ip_error = np.full(len(self._times), np.nan)
        dip_dt = np.full(len(self._times), np.nan)
        dipprog_dt = np.full(len(self._times), np.nan)
        # Get measured plasma current parameters
        try:
            t_ip = self.conn.get(
                f"dim_of(ptdata('ip', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip = self.conn.get(f"ptdata('ip', {self._shot_id})").data()  # [A]
            dip_dt = np.gradient(ip, t_ip)
            ip = interp1(t_ip, ip, self._times, 'linear')
            dip_dt = interp1(t_ip, dip_dt, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get measured plasma current parameters")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Get programmed plasma current parameters
        try:
            t_ip_prog = self.conn.get(
                f"dim_of(ptdata('iptipp', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_prog = self.conn.get(
                f"ptdata('iptipp', {self._shot_id})").data()  # [A]
            polarity = np.unique(self.conn.get(
                f"ptdata('iptdirect', {self._shot_id})").data())
            if len(polarity) > 1:
                self.logger.info(
                    f"[Shot {self._shot_id}]:Polarity of Ip target is not constant. Using value at first timestep.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]: Polarity array {polarity}")
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, self._times, 'linear')
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get programmed plasma current parameters")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Now get the signal pointname 'ipimode'.  This PCS signal denotes whether
        # or not PCS is actually feedback controlling the plasma current.  There
        # are times when feedback of Ip is purposely turned off, such as during
        # electron cyclotron current drive experiments.  Here is how to interpret
        # the value of 'ipimode':
        #  0: normal Ip feedback to E-coils supplies
        #  3: almost normal Ip feedback, except that abs(Ip) > 2.5 MA
        #  Anything else: not in normal Ip feedback mode.  In this case, the
        # 'ip_prog' signal is irrelevant, and therefore 'ip_error' is not defined.
        try:
            ipimode = self.conn.get(
                f"ptdata('ipimode', {self._shot_id})").data()
            t_ipimode = self.conn.get(
                f"dim_of(ptdata('ipimode', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ipimode = interp1(t_ipimode, ipimode, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get ipimode signal. Setting to NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            ipimode = np.full(len(self._times), np.nan)
        feedback_on_indices = np.where((ipimode == 0) | (ipimode == 3))
        ip_error[feedback_on_indices] = ip[feedback_on_indices] - \
            ip_prog[feedback_on_indices]
        # Finally, get 'epsoff' to determine if/when the E-coil power supplies have railed
        # Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
        # PCS feedback control of Ip is not being applied.  Therefore the
        # 'ip_error' parameter is undefined for these times.
        try:
            epsoff = self.conn.get(f"ptdata('epsoff', {self._shot_id})").data()
            t_epsoff = self.conn.get(
                f"dim_of(ptdata('epsoff', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            t_epsoff += .001  # Avoid problem with simultaneity of epsoff being triggered exactly on the last time sample
            epsoff = interp1(t_epsoff, epsoff, self._times, 'linear')
            railed_indices = np.where(np.abs(epsoff) > .5)
            power_supply_railed = np.zeros(len(self._times))
            power_supply_railed[railed_indices] = 1
            ip_error[railed_indices] = np.nan
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get epsoff signal. Setting to NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            power_supply_railed = np.full(len(self._times), np.nan)
        # 'ip_prog': ip_prog,
        return pd.DataFrame({'ip': ip, 'ip_error': ip_error, 'dip_dt': dip_dt, 'dipprog_dt': dipprog_dt, 'power_supply_railed': power_supply_railed})

    def get_rt_ip_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        ip_rt = np.full(len(self._times), np.nan)
        ip_prog_rt = np.full(len(self._times), np.nan)
        ip_error_rt = np.full(len(self._times), np.nan)
        dip_dt_rt = np.full(len(self._times), np.nan)
        dipprog_dt_rt = np.full(len(self._times), np.nan)
        # Get measured plasma current parameters
        try:
            t_ip_rt = self.conn.get(
                f"dim_of(ptdata('ipsip', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_rt = self.conn.get(
                f"ptdata('ipsip', {self._shot_id})").data()*1.e6  # [MA] -> [A]
            dip_dt_rt = np.gradient(ip_rt, t_ip_rt)
            ip_rt = interp1(t_ip_rt, ip_rt, self._times, 'linear')
            dip_dt_rt = interp1(t_ip_rt, dip_dt_rt, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get measured plasma current parameters")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Get programmed plasma current parameters
        try:
            t_ip_prog_rt = self.conn.get(
                f"dim_of(ptdata('ipsiptargt', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_prog_rt = self.conn.get(
                f"ptdata('ipsiptargt', {self._shot_id})").data()*1.e6*.5  # [MA] -> [A]
            polarity = np.unique(self.conn.get(
                f"ptdata('iptdirect', {self._shot_id})").data())
            if len(polarity) > 1:
                self.logger.info(
                    f"[Shot {self._shot_id}]:Polarity of Ip target is not constant. Setting to first value in array.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]: Polarity array: {polarity}")
                polarity = polarity[0]
            ip_prog_rt = ip_prog_rt * polarity
            dipprog_dt_rt = np.gradient(ip_prog_rt, t_ip_prog_rt)
            ip_prog_rt = interp1(t_ip_prog_rt, ip_prog_rt,
                                 self._times, 'linear')
            dipprog_dt_rt = interp1(
                t_ip_prog_rt, dipprog_dt_rt, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get programmed plasma current parameters")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        try:
            t_ip_error_rt = self.conn.get(
                f"dim_of(ptdata('ipeecoil', {self._shot_id}))").data()/1.e3  # [ms] to [s]
            ip_error_rt = self.conn.get(
                f"ptdata('ipeecoil', {self._shot_id})").data()*1.e6*.5  # [MA] -> [A]
            ip_error_rt = interp1(
                t_ip_error_rt, ip_error_rt, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get ipeecoil signal. Setting to NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Now get the signal pointname 'ipimode'.  This PCS signal denotes whether
        # or not PCS is actually feedback controlling the plasma current.  There
        # are times when feedback of Ip is purposely turned off, such as during
        # electron cyclotron current drive experiments.  Here is how to interpret
        # the value of 'ipimode':
        #  0: normal Ip feedback to E-coils supplies
        #  3: almost normal Ip feedback, except that abs(Ip) > 2.5 MA
        #  Anything else: not in normal Ip feedback mode.  In this case, the
        # 'ip_prog' signal is irrelevant, and therefore 'ip_error' is not defined.
        try:
            ipimode = self.conn.get(
                f"ptdata('ipimode', {self._shot_id})").data()
            t_ipimode = self.conn.get(
                f"dim_of(ptdata('ipimode', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ipimode = interp1(t_ipimode, ipimode, self._times, 'linear')
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get ipimode signal. Setting to NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            ipimode = np.full(len(self._times), np.nan)
        feedback_off_indices = np.where((ipimode != 0) & (ipimode == 3))
        ip_error_rt[feedback_off_indices] = np.nan
        # Finally, get 'epsoff' to determine if/when the E-coil power supplies have railed
        # Times at which power_supply_railed ~=0 (i.e. epsoff ~=0) mean that
        # PCS feedback control of Ip is not being applied.  Therefore the
        # 'ip_error' parameter is undefined for these times.
        try:
            epsoff = self.conn.get(f"ptdata('epsoff', {self._shot_id})").data()
            t_epsoff = self.conn.get(
                f"dim_of(ptdata('epsoff', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            t_epsoff += .001  # Avoid problem with simultaneity of epsoff being triggered exactly on the last time sample
            epsoff = interp1(t_epsoff, epsoff, self._times, 'linear')
            railed_indices = np.where(np.abs(epsoff) > .5)
            power_supply_railed = np.zeros(len(self._times))
            power_supply_railed[railed_indices] = 1
            ip_error_rt[railed_indices] = np.nan
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get epsoff signal. power_supply_railed will be NaN.")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            power_supply_railed = np.full(len(self._times), np.nan)
        # 'dip_dt_RT': dip_dt_rt,
        return pd.DataFrame({'ip_RT': ip_rt, 'ip_error_RT': ip_error_rt, 'dipprog_dt_RT': dipprog_dt_rt, 'power_supply_railed': power_supply_railed})

    def get_z_parameters(self):
        """
        On DIII-D the plasma control system uses isoflux
        control to control the plasma shape and position.  It does
        NOT use zcur control.  Therefore, the PCS does not have a
        programmed vertical position.  This this routine will now
        always return an arrays of NaN for z_prog, z_error, and
        z_error_norm.
        """
        z_cur = np.full(len(self._times), np.nan)
        z_cur_norm = np.full(len(self._times), np.nan)
        z_prog = np.full(len(self._times), np.nan)
        z_error = np.full(len(self._times), np.nan)
        z_error_norm = np.full(len(self._times), np.nan)
        self.conn.openTree('d3d', self._shot_id)
        try:
            t_z_cur = self.conn.get(
                f"dim_of(ptdata('vpszp', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            z_cur = self.conn.get(
                f"ptdata('vpszp', {self._shot_id})").data()/1.e2  # [cm] -> [m]
            z_cur = interp1(t_z_cur, z_cur, self._times, 'linear')
            self.conn.openTree(self.efit_tree_name, self._shot_id)
            try:
                t_a = self.conn.get(
                    r'\efit_a_eqdsk:atime').data()/1.e3  # [ms] -> [s]
                a_minor = self.conn.get(r'\efit_a_eqdsk:aminor').data()  # [m]
                chisq = self.conn.get(r'\efit_a_eqdsk:chisq').data()
                invalid_indices = np.where(chisq > 50)
                a_minor[invalid_indices] = np.nan
                a_minor = interp1(t_a, a_minor, self._times, 'linear')
                z_cur_norm = z_cur/a_minor
            except MdsException as e:
                self.logger.info(
                    f"[Shot {self._shot_id}]:Failed to get efit parameters")
                self.logger.debug(
                    f"[Shot {self._shot_id}]:{traceback.format_exc()}")
                z_cur_norm = z_cur / self.nominal_flattop_radius
        except MdsException as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get vpszp signal")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        return pd.DataFrame({'zcur': z_cur, 'zcur_normalized': z_cur_norm, 'z_prog': z_prog, 'z_error': z_error, 'z_error_normalized': z_error_norm})

    def get_n1_bradial_parameters(self):
        # The following shots are missing bradial calculations in MDSplus and must be loaded from a separate datafile
        if self._shot_id >= 176030 and self._shot_id <= 176912:
            # TODO: Move to a folder like "/fusion/projects/disruption_warning/data"
            filename = '/fusion/projects/disruption_warning/matlab_programs/recalc.nc'
            ncid = nc.Dataset(filename, 'r')
            brad = ncid.variables['dusbradial_calculated'][:]
            t_n1 = ncid.variables['times'][:]*1.e-3  # [ms] -> [s]
            shots = ncid.variables['shots'][:]
            shot_indices = np.where(shots == self._shot_id)
            if len(shot_indices) == 1:
                dusbradial = brad[shot_indices, :]*1.e-4  # [T]
            else:
                self.logger.info(
                    f"Shot {self._shot_id} not found in {filename}.  Returning NaN.")
                dusbradial = np.full(len(self._times), np.nan)
            ncid.close()
        # Check ONFR than DUD(legacy)
        else:
            try:
                dusbradial, t_n1 = self._get_signal(
                    f"ptdata('onsbradial',{self._shot_id})")*1.e-4  # [T]
            except MdsException as e:
                self.logger.debug(
                    f"[Shot {self._shot_id}]:{traceback.format_exc()}")
                try:
                    dusbradial, t_n1 = self._get_signal(
                        f"ptdata('dusbradial',{self._shot_id})")*1.e-4  # [T]
                except MdsException as e:
                    self.logger.info(
                        f"[Shot {self._shot_id}]:Failed to get n1 bradial signal. Returning NaN.")
                    self.logger.debug(
                        f"[Shot {self._shot_id}]:{traceback.format_exc()}")
                    n_equal_1_mode = np.full(len(self._times), np.nan)
                    n_equal_1_normalized = np.full(len(self._times), np.nan)
                    return pd.DataFrame({'n_equal_1_normalized': n_equal_1_normalized, 'n_equal_1_mode': n_equal_1_mode})
        n_equal_1_mode = interp1(dusbradial, t_n1, self._times)
        # Get toroidal field Btor
        b_tor, _ = self._get_signal(
            "ptdata('bt',{self._shot_id})")  # [T]
        n_equal_1_normalized = n_equal_1_mode/b_tor
        return pd.DataFrame({'n_equal_1_normalized': n_equal_1_normalized, 'n_equal_1_mode': n_equal_1_mode})

    def get_n1rms_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        n1rms, t_n1rms = self._get_signal(r'\n1rms', interpolate=False)
        n1rms *= 1.e-4  # Gauss -> Tesla
        n1rms = interp1(t_n1rms, n1rms, self._times)
        b_tor = self._get_signal(
            "ptdata('bt',{self._shot_id})")  # [T]
        n1rms_norm = n1rms / np.abs(b_tor)
        return pd.DataFrame({'n1rms': n1rms, 'n1rms_normalized': n1rms_norm})

    # TODO: Need to test and unblock recalculating peaking factors
    # By default get_peaking_factors should grab the data from MDSPlus as opposed to recalculate. See DPP v4 document for details:
    # https://docs.google.com/document/d/1R7fI7mCOkMQGt8xX2nS6ZmNNkcyvPQ7NmBfRPICFaFs/edit?usp=sharing
    def get_peaking_factors(self):
        ts_data_type = 'blessed'  # either 'blessed', 'unblessed', or 'ptdata'
        # metric to use for core/edge binning (either 'psin' or 'rhovn')
        ts_radius = 'rhovn'
        # ts_radius value defining boundary of 'core' region (between 0 and 1)
        ts_core_margin = 0.3
        # All data outside this range excluded. For example, psin=0 at magnetic axis and 1 at separatrix.
        ts_radial_range = (0, 1)
        # set to true to interpolate ts_channel data onto equispaced radial grid
        ts_equispaced = False
        # fan to use for P_rad peaking factors (either 'lower', 'upper', or 'custom')
        bolometer_fan = 'custom'
        # array of bolometer fan channel numbers covering divertor (upper fan: 1->24, lower fan: 25:48)
        div_channels = np.arange(3, 8)+24
        # time window for filtering raw bolometer signal in [ms]
        smoothing_window = 40
        p_rad_core_def = 0.06  # percentage of DIII-D veritcal extent defining the core margin
        # 'brightness'; % either 'brightness' or 'power' ('z')
        p_rad_metric = 'brightness'
        # Ts options
        ts_options = ['combined', 'core', 'tangential']
        # vertical range of the DIII-D cross section in meters
        vert_range = 3
        te_pf = np.full(len(self._times), np.nan)
        ne_pf = np.full(len(self._times), np.nan)
        rad_cva = np.full(len(self._times), np.nan)
        rad_xdiv = np.full(len(self._times), np.nan)
        try:
            rad_cva = self._get_signal(
                f"ptdata('dpsradcva', {self._shot_id})")
            rad_xdiv = self._get_signal(
                f"ptdata('dpsradxdiv', {self._shot_id})")
        except MdsException as e:
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get CVA and XDIV from MDSPlus. Calculating locally, results may be inaccurate.")
            rad_cva = np.full(len(self._times), np.nan)
            rad_xdiv = np.full(len(self._times), np.nan)
        try:
            ts = self._get_ne_te()
            for option in ts_options:
                if option in ts:
                    ts = ts[option]
            efit_dict = self._get_efit_dict()
        except Exception as e:
            self.logger.info(f"[Shot {self._shot_id}]:Failed to get TS data")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            ts = 0
        try:
            ts['psin'], ts['rhovn'] = efit_rz_interp(ts, efit_dict)
            print(ts['rhovn'].shape)
        except Exception as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to interpolate TS data")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        try:
            p_rad = self._get_p_rad()
        except Exception as e:
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to get bolometer data")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            p_rad = 0
        if p_rad == 0 and ts == 0:
            self.logger.info(
                f"[Shot {self._shot_id}]:Both TS and bolometer data missing for shot")
        if ts != 0 and ts_radius in ts:
            # Drop data outside of valid range
            invalid_indices = np.where((ts[ts_radius] < ts_radial_range[0]) | (
                ts[ts_radius] > ts_radial_range[1]))
            print(ts['te'].shape)
            ts['te'][invalid_indices] = np.nan
            ts['ne'][invalid_indices] = np.nan
            ts['te'][np.isnan(ts[ts_radius])] = np.nan
            ts['ne'][np.isnan(ts[ts_radius])] = np.nan
            if ts_equispaced:
                raise NotImplementedError(
                    "Equispaced is currently assumed to be false")  # TODO
            # Find core bin for Thomson and calculate Te, ne peaking factors
            core_mask = ts[ts_radius] < ts_core_margin
            te_core = ts['te']
            te_core[~core_mask] = np.nan
            ne_core = ts['ne']
            ne_core[~core_mask] = np.nan
            te_pf = np.nanmean(te_core, axis=1)/np.nanmean(ts['te'], axis=1)
            ne_pf = np.nanmean(ne_core, axis=1)/np.nanmean(ts['ne'], axis=1)
            # Calculate Prad CVA, X-DIV Peaking Factors
            # # Interpolate zmaxis and channel intersects x onto the bolometer timebase
            z_m_axis = interp1(efit_dict['time'],
                               efit_dict['zmaxis'], ts['time'])
            z_m_axis = np.repeat(
                z_m_axis[:, np.newaxis], p_rad['x'].shape[1], axis=1)
            p_rad['xinterp'] = interp1(p_rad['xtime'], p_rad['x'], p_rad['t'])
            # # Determine the bolometer channels falling in the 'core' bin
            core_indices = (p_rad['xinterp'] < z_m_axis + p_rad_core_def*vert_range) & (
                p_rad['xinterp'] > z_m_axis - p_rad_core_def*vert_range)
            # # Designate the divertor bin and find all 'other' channels not in that bin
            div_indices = np.searchsorted(p_rad['ch_avail'], div_channels)
            other_indices = ~div_indices
            # # Grab p_rad measurements for each needed set of channels
            p_rad_core = p_rad[p_rad_metric]
            p_rad_all_but_core = p_rad_core.copy()
            # QUESTION: Why fill with nans for core but just keep valid indices for divertor
            p_rad_core[~core_indices] = np.nan
            p_rad_all_but_core[core_indices] = np.nan
            p_rad_div = p_rad[p_rad_metric][div_indices, :]
            p_rad_all_but_div = p_rad[p_rad_metric][other_indices, :]
            # # Calculate the peaking factors
            rad_cva = np.nanmean(p_rad_core, axis=0) / \
                np.nanmean(p_rad_all_but_div, axis=0)
            rad_xdiv = np.nanmean(p_rad_div, axis=0) / \
                np.nanmean(p_rad_all_but_core, axis=0)
            rad_cva = interp1(p_rad['t'], rad_cva, self._times)
            rad_xdiv = interp1(p_rad['t'], rad_xdiv, self._times)
        return pd.DataFrame({'te_pf': te_pf, 'ne_pf': ne_pf, 'rad_cva': rad_cva, 'rad_xdiv': rad_xdiv})

    # TODO: Finish implementing just in case
    def _efit_map_rz_to_rho_original(self, ts_dict, efit_dict):
        slices = np.zeros(ts_dict['time'].shape)
        # If thomson starts before EFIT (often does), then use the first valid EFIT slice for early Thomson data.
        early_indices = np.where(ts_dict['time'] < efit_dict['time'])
        if len(early_indices[0]) > 0:
            slices[early_indices] = 1
            first_ts = early_indices[0][-1]
        else:
            first_ts = 0
        # If Thomson ends after EFIT (also often happens), then use the last valid EFIT slice for late Thomson data.
        late_indices = np.where(ts_dict['time'] >= efit_dict['time'])
        if len(late_indices[0]) > 0:
            slices[late_indices] = len(efit_dict['time'])
            last_ts = late_indices[0][0] - 1
        else:
            last_ts = len(ts_dict['time']) - 1
        diag_slices = np.arange(first_ts, last_ts+1, 1)
        # Acquire list of diag time slices w/in EFIT time range; Should find closest EFIT for each one
        for i in diag_slices:
            slices[i] = np.argmin(
                np.abs(efit_dict['time'] - ts_dict['time'][i]))
        # Interpolate EFIT data onto Thomson time slices
        psin_diag_arr = np.zeros((len(efit_dict['time']), len(ts_dict['z'])))
        for r in np.unique(ts_dict['r']):
            dr = r - efit_dict['r']
            # Find closet EFIT R on the left and right
            right = np.where(efit_dict['r'] > r, 1)
            left = right - 1
            if efit_dict['r'][right] == r:
                psin_slice = np.squeeze(efit_dict['psin'][:, right, :])

    def get_core_edge_vals(self):
        ##################################################
        # Settings
        ts_data_type = 'blessed'  # either 'blessed', 'unblessed', or 'ptdata'
        # metric to use for core/edge binning (either 'psin' or 'rhovn')
        ts_radius = 'rhovn'
        # ts_radius value defining boundary of 'core' region (between 0 and 1)
        ts_core_margin = 0.3
        # ts_radius value defining inner and outer side of 'edge' region (between ts_core_margin and 1)
        ts_edge_inner = 0.85
        ts_edge_outer = 0.95
        # All data outside this range excluded. For example, psin=0 at magnetic axis and 1 at separatrix.
        ts_radial_range = (0, 1)
        # set to true to interpolate ts_channel data onto equispaced radial grid
        ts_equispaced = True
        ###################################################

        # Initialize arrays
        te_core = np.full(len(self._times), np.nan)
        ne_core = np.full(len(self._times), np.nan)
        # Averaged over edge region
        te_edge = np.full(len(self._times), np.nan)
        ne_edge = np.full(len(self._times), np.nan)
        # Averaged over 85th to 88th surface
        te_edge_80to85 = np.full(len(self._times), np.nan)
        ne_edge_80to85 = np.full(len(self._times), np.nan)
        te_edge_85to90 = np.full(len(self._times), np.nan)
        ne_edge_85to90 = np.full(len(self._times), np.nan)
        te_edge_90to95 = np.full(len(self._times), np.nan)
        ne_edge_90to95 = np.full(len(self._times), np.nan)
        te_edge_95to100 = np.full(len(self._times), np.nan)
        ne_edge_95to100 = np.full(len(self._times), np.nan)
        # Separatrix
        te_sep = np.full(len(self._times), np.nan)
        ne_sep = np.full(len(self._times), np.nan)

        # Try to get data via _get_ne_te()
        try:
            ts = self._get_ne_te()
            efit_dict = self._get_efit_dict
            ts['psin'], ts['rhovn'] = efit_rz_interp(ts, efit_dict)
        except Exception as e:
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            ts = 0
        if ts == 0:
            self.logger.info(
                f"[Shot {self._shot_id}]:Both TS data missing for shot #{self._shot_id}")
        if ts != 0:
            # Drop data outside of valid range #ADM: this looks unfinished
            invalid_indices = np.where((ts[ts_radius] < ts_radial_range[0]) | (
                ts[ts_radius] > ts_radial_range[1]))

        # TODO: 1) Interpolate in core and edge regions, 2) compute average in these regions and store in respective array. Note that we may need to expand the available indices beyond 1

        return pd.DataFrame({'te_core': te_core, 'ne_core': ne_core, 'te_core': te_edge, 'ne_edge': ne_edge, 'te_edge_80to85': te_edge_80to85, 'ne_edge_80to85': ne_edge_80to85,
                             'te_edge_85to90': te_edge_85to90, 'ne_edge_85to90': ne_edge_85to90, 'te_edge_90to95': te_edge_90to95, 'ne_edge_90to95': ne_edge_90to95, 'te_edge_95to100': te_edge_95to100, 'ne_edge_95to100': ne_edge_95to100, 'te_sep': te_sep, 'ne_sep': ne_sep})

    def get_zeff_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        # Get Zeff
        try:
            zeff = self.conn.get(
                r'\d3d::top.spectroscopy.vb.zeff:zeff').data()
            # t_nbi = self.conn.get(
            # r"dim_of(\d3d::top.nb:pinj)").data()/1.e3  # [ms]->[s]
            t_zeff = self.conn.get(
                r'dim_of(\d3d::top.spectroscopy.vb.zeff:zeff)').data()/1.e3  # [ms] -> [s]
            if len(t_zeff) > 2:
                zeff = interp1(t_zeff, zeff, self._times,
                               'linear', bounds_error=False, fill_value=0.)
            else:
                zeff = np.zeros(len(self._times))
                self.logger.info(
                    f"[Shot {self._shot_id}]:No zeff data found in this shot.")
        except MdsException as e:
            zeff = np.zeros(len(self._times))
            self.logger.info(
                f"[Shot {self._shot_id}]:Failed to open Zeff node")
            self.logger.debug(
                f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        return pd.DataFrame({'z_eff': zeff})

    def get_kappa_area(self):
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        a_minor = self.conn.get(r'\efit_a_eqdsk:aminor').data()
        area = self.conn.get(r'\efit_a_eqdsk:area').data()
        chisq = self.conn.get(r'\efit_a_eqdsk:chisq').data()
        t = self.conn.get(r'\efit_a_eqdsk:atime')
        kappa_area = area / (np.pi * a_minor**2)
        invalid_indices = np.where(chisq > 50)
        kappa_area[invalid_indices] = np.nan
        kappa_area = interp1(t, kappa_area, self._times)
        return pd.DataFrame({'kappa_area': kappa_area})

    def get_h_parameters(self):
        h98 = np.full(len(self._times), np.nan)
        self.conn.openTree('transport', self._shot_id)
        h98, t_h98 = self._get_signal(r'\H_THH98Y2')
        self.conn.openTree('d3d')
        h_alpha, t_h_alpha = self._get_signal(r'\fs04')
        h98 = interp1(t_h98, h98, self._times)
        h_alpha = interp1(t_h_alpha, h_alpha, self._times)
        return pd.DataFrame({'H98': h98, 'H_alpha': h_alpha})

    def get_shape_parameters(self):
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        efit_time = self.conn.get(
            r'\efit_a_eqdsk:atime').data()/1.e3  # [ms] -> [s]
        sqfod = self.conn.get(r'\efit_a_eqdsk:sqfod').data()
        sqfou = self.conn.get(r'\efit_a_eqdsk:sqfou').data()
        tritop = self.conn.get(r'\efit_a_eqdsk:tritop').data()  # meters
        tribot = self.conn.get(r'\efit_a_eqdsk:tribot').data()  # meters
        # plasma minor radius [m]
        aminor = self.conn.get(r'\efit_a_eqdsk:aminor').data()
        chisq = self.conn.get(r'\efit_a_eqdsk:chisq').data()
        # Compute triangularity and squareness:
        delta = (tritop+tribot)/2.0
        squareness = (sqfod+sqfou)/2.0

        # Remove invalid indices
        invalid_indices = np.where(chisq > 50)
        delta[invalid_indices] = np.nan
        squareness[invalid_indices] = np.nan
        aminor[invalid_indices] = np.nan

        # Interpolate to desired times
        delta = interp1(efit_time, delta, self._times, 'linear',
                        bounds_error=False, fill_value=np.nan)
        squareness = interp1(efit_time, squareness, self._times,
                             'linear', bounds_error=False, fill_value=np.nan)
        aminor = interp1(efit_time, aminor, self._times,
                         'linear', bounds_error=False, fill_value=np.nan)
        return pd.DataFrame({'delta': delta, 'squareness': squareness, 'aminor': aminor})

    def _get_ne_te(self, data_source="blessed", ts_systems=['core', 'tangential']):
        if data_source == 'blessed':  # 'blessed' by Thomson group
            mds_path = r'\top.ts.blessed.'
        elif data_source == 'unblessed':
            mds_path = r'\top.ts.revisions.revision00.'
        elif data_source == 'ptdata':
            mds_path = r'\top.ts.blessed.'  # Don't ask...I don't have the answer
            raise NotImplementedError(
                "ptdata case not fully implemented yet")  # TODO
        else:
            raise ValueError(f"Invalid data_source: {data_source}")
        # Account for pointname formatting change in 2017 (however using ptdata is unimplemented)
        suffix = {'core': 'cor', 'tangential': 'tan'}
        if self._shot_id < 172749:  # First shot on Sep 19, 2017
            suffix['tangential'] = 'hor'
        self.conn.openTree('electrons', self._shot_id)
        lasers = dict()
        for laser in ts_systems:
            lasers[laser] = dict()
            sub_tree = f"{mds_path}{laser}"
            try:
                lasers[laser]['time'] = self.conn.get(
                    f"dim_of({sub_tree}:temp,0)").data()/1.e3  # [ms] -> [s]
            except MdsException as e:
                lasers[laser] = None
                self.logger.info(
                    f"[Shot {self._shot_id}]: Failed to get {laser} time. Setting laser data to None.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]:{traceback.format_exc()}")
                continue
            child_nodes = {'r': 'r', 'z': 'z', 'te': 'temp', 'ne': 'density',
                           'time': 'time', 'te_error': 'temp_e', 'ne_error': 'density_e'}
            for node, name in child_nodes.items():
                try:
                    lasers[laser][node] = self.conn.get(
                        f"{sub_tree}:{name}").data()
                except MdsException as e:
                    lasers[laser][node] = np.full(
                        lasers[laser]['time'].shape, np.nan)
                    self.logger.info(
                        f"[Shot {self._shot_id}]: Failed to get {laser}:{name}({node}) data, Setting to all NaNs.")
                    self.logger.debug(
                        f"[Shot {self._shot_id}]:{traceback.format_exc()}")
            # Place NaNs for broken channels
            lasers[laser]['te'][np.where(
                lasers[laser]['te'] == 0)] = np.nan
            lasers[laser]['ne'][np.where(
                lasers[laser]['ne'] == 0)] = np.nan
        # If both systems/lasers available, combine them and interpolate the data
        # from the tangential system onto the finer (core) timebase
        if 'tangential' in lasers and lasers['tangential'] is not None:
            if 'core' in lasers and lasers['core'] is not None:
                lasers['combined'] = dict()
                # Interpolate tangential data onto core timebase
                for key in lasers['tangential']:
                    if key not in ['time', 'r', 'z']:
                        lasers['tangential'][key] = interp1(
                            lasers['tangential']['time'], lasers['tangential'][key], lasers['core']['time'])
                        lasers['combined'][key] = np.concatenate((
                            lasers['core'][key], lasers['tangential'][key]))
                lasers['tangential']['time'] = lasers['core']['time']
                lasers['combined']['time'] = lasers['core']['time']
                lasers['combined']['r'] = np.concatenate(
                    (lasers['core']['r'], lasers['tangential']['r']))
                lasers['combined']['z'] = np.concatenate(
                    (lasers['core']['z'], lasers['tangential']['z']))
        return lasers

    def _get_p_rad(self, fan='custom'):
        if fan == 'upper':
            fan_chans = np.arange(0, 24)
        elif fan == 'lower':
            fan_chans = np.arange(24, 48)
        elif fan == 'custom':
            # 1st choice (heavily cover divertor and core)
            fan_chans = np.array(
                [3, 4, 5, 6, 7, 8, 9, 12, 14, 15, 16, 22]) + 24

        # Get bolometry data
        self.conn.openTree("bolom", self._shot_id)
        bol_prm, _ = self._get_signal(r'\bol_prm', interpolate=False)
        lower_channels = [f"bol_u{i+1:02d}_v" for i in range(24)]
        upper_channels = [f"bol_l{i+1:02d}_v" for i in range(24)]
        bol_channels = lower_channels + upper_channels
        bol_signals = []
        bol_times = []  # TODO: Decide whether to actually use all bol_times instead of just first one
        for i in range(48):
            bol_signal, bol_time = self._get_signal(
                fr"\top.raw:{bol_channels[i]}", interpolate=False)
            bol_signals.append(bol_signal)
            bol_times.append(bol_time)
        a_struct = get_bolo(self._shot_id, bol_channels,
                            bol_prm, bol_signals, bol_times[0])
        b_struct = power(a_struct)
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        r_major_axis, efit_time = self._get_signal(
            r'\top.results.geqdsk:rmaxis', interpolate=False)
        data_dict = {'ch_avail': [], 'z': [], 'brightness': [],
                     'power': [], 'x': np.full((len(efit_time), len(fan_chans)), np.nan), 'xtime': efit_time, 't': a_struct.raw_time}
        for i in range(len(fan_chans)):
            chan = fan_chans[i]
            data_dict['power'].append(b_struct.chan[chan].chanpwr)
            if a_struct.channels[chan].ier == 0:
                data_dict['ch_avail'].append(chan)
            data_dict['x'][:, i] = a_struct.channels[chan].Z + \
                np.tan(a_struct.channels[chan].angle*np.pi/180.0) * \
                (r_major_axis - a_struct.channels[chan].R)
            b_struct.chan[chan].chanpwr[np.where(
                b_struct.chan[chan].chanpwr < 0)] = 0
            b_struct.chan[chan].brightness[np.where(
                b_struct.chan[chan].brightness < 0)] = 0
            data_dict['z'].append(b_struct.chan[i].chanpwr)
            data_dict['brightness'].append(b_struct.chan[i].brightness)
        return data_dict

    # TODO: Replace all instances of efit_dict with a dataclass
    def _get_efit_dict(self):
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        efit_dict = dict()
        path = r'\top.results.geqdsk:'
        nodes = ['z', 'r', 'rhovn', 'psirz', 'zmaxis', 'ssimag', 'ssibry']
        efit_dict['time'] = self.conn.get(
            f"dim_of({path}psirz,2)").data()/1.e3  # [ms] -> [s]
        for node in nodes:
            try:
                efit_dict[node] = self.conn.get(f"{path}{node}").data()
            except MdsException as e:
                efit_dict[node] = np.full(efit_dict['time'].shape, np.nan)
                self.logger.info(
                    f"[Shot {self._shot_id}]: Failed to get {node} from efit, Setting to all NaNs.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        # Normalize the poloidal flux grid (0=magnetic axis, 1=boundary)
        # [Translated from D. Eldon's OMFITeqdsk read_basic_eq_from_mds() function]
        psi_norm_f = efit_dict['ssibry'] - efit_dict['ssimag']
        problems = np.where(psi_norm_f == 0)[0]
        # Prevent divide by 0 error by replacing 0s in the denominator
        psi_norm_f[problems] = 1
        efit_dict['psin'] = (efit_dict['psirz'] - efit_dict['ssimag']
                             [:, np.newaxis, np.newaxis])/psi_norm_f[:, np.newaxis, np.newaxis]
        efit_dict['psin'][problems, :, :] = 0
        return efit_dict


if __name__ == '__main__':
    logger = logging.getLogger('disruption_py')
    logger.setLevel(logging.DEBUG)
    shot = D3DShot(D3D_DISRUPTED_SHOT, 'EFIT05',
                   disruption_time=4.369214483261109, populate='l_mode')
    print(shot.data.columns)
    # print(shot.data[['te_pf','ne_pf','rad_cva','rad_xdiv']])
    print(shot.data.head())
