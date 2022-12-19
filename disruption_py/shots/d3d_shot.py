from shot import Shot
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

import pandas as pd
import numpy as np

import MDSplus
from MDSplus import *

from disruption_py.utils import interp1, interp2, smooth, gaussian_fit, gsastd, get_bolo, power
import disruption_py.data
D3D_DISRUPTED_SHOT = 175552
"""
Useful Examples:
https://diii-d.gat.com/diii-d/MDSplusAPI_Python_pg1
https://diii-d.gat.com/diii-d/Gadata_py
"""


class D3DShot(Shot):
    efit_vars = {'beta_n': '\efit_a_eqdsk:betan', 'beta_p': '\efit_a_eqdsk:betap', 'kappa': '\efit_a_eqdsk:kappa', 'li': '\efit_a_eqdsk:li', 'upper_gap': '\efit_a_eqdsk:gaptop', 'lower_gap': '\efit_a_eqdsk:gapbot',
                 'q0': '\efit_a_eqdsk:q0', 'qstar': '\efit_a_eqdsk:qstar', 'q95': '\efit_a_eqdsk:q95',  'v_loop_efit': '\efit_a_eqdsk:vsurf', 'Wmhd': '\efit_a_eqdsk:wmhd', 'chisq': '\efit_a_eqdsk:chisq'}
    efit_derivs = ['beta_p', 'li', 'Wmhd']
    nominal_flattop_radius = 0.59

    def __init__(self, shot_id, efit_tree_name, data_columns=DEFAULT_SHOT_COLUMNS, data=None):
        super().__init__(shot_id, data_columns, data)
        self.conn = MDSplus.Connection('atlas.gat.com')
        self.efit_tree_name = str(efit_tree_name)
        if data is None:
            self.data = pd.DataFrame()
            self._populate_shot_data()
        self._times = None  # TODO: Set somehow

    def _populate_shot_data(self):
        efit_data = self.get_efit_parameters()
        self.data = pd.concat([self.data, efit_data], ignore_index=True)

    def get_efit_parameters(self):
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        efit_data = {k: self.conn.get(v).data()
                     for k, v in self.efit_vars.items()}
        efit_time = self.conn.get('\efit_a_eqdsk:atime').data()/1000.0
        self._times = efit_time  # TODO: Reconsider how shot times are chosen
        # EFIT reconstructions are sometimes invalid, particularly when very close
        # to a disruption.  There are a number of EFIT parameters that can indicate
        # invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
        # 'chisq' to determine which time slices should be excluded from our
        # disruption warning database.
        invalid_indices = np.where(efit_data['chisq'] > 50)
        for param in efit_data:
            efit_data[param][invalid_indices] = np.nan
        self._times[invalid_indices] = np.nan
        for param in self.efit_derivs:
            efit_data['d' + param +
                      '_dt'] = np.gradient(efit_data[param], efit_time)
        if not np.array_equal(self._times, efit_time):
            for param in efit_data:
                efit_data[param] = interp1(
                    efit_time, efit_data[param], self._times)
        return pd.DataFrame([efit_data])

    def get_power_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        # Get neutral beam injected power
        try:
            p_nbi = self.conn.get(
                r"\d3d::top.nb:pinj").data()*1.e3  # [KW]->[W]
            t_nbi = self.conn.get(
                r"dim_of(\d3d::top.nb:pinj)").data()*1.e3  # [ms]->[s]
            if len(t_nbi) > 2:
                p_nbi = interp1(t_nbi, p_nbi, self._times,
                                'linear', bounds_error=False, fill_value=0.)
            else:
                print("No NBI power data found in this shot.")
                p_nbi = np.zeros(len(self._times))
        except mdsExceptions.TreeFOPENR as e:
            p_nbi = np.zeros(len(self._times))
            print("Failed to open NBI node")
        # Get electron cycholotrn heating (ECH) power. It's poitn data, so it's not stored in an MDSplus tree
        self.conn.openTree('rf', self._shot_id)
        try:
            p_ech = self.conn.get(r"\top.ech.totl:echpwrc").data()  # [W]
            t_ech = self.conn.get(
                r"dim_of(\top.ech.totl:echpwrc)").data()*1.e3  # [ms]->[s]
            if len(t_ech) > 2:
                p_ech = interp1(t_ech, p_ech, self._times,
                                'linear', bounds_error=False, fill_value=0.)
            else:
                print("No ECH power data found in this shot.")
                p_ech = np.zeros(len(self._times))
        except mdsExceptions.TreeFOPENR as e:
            p_ech = np.zeros(len(self._times))
            print("Failed to open ECH node")
        # Get ohmic power and loop voltage
        p_ohm, v_loop = self.get_ohmic_parameters()
        # Radiated power
        # We had planned to use the standard signal '\bolom::prad_tot' for this
        # parameter.  However, the processing involved in calculating \prad_tot
        # from the arrays of bolometry channels involves non-causal filtering with
        # a 50 ms window.  This is not acceptable for our purposes.  Tony Leonard
        # provided us with the two IDL routines that are used to do the automatic
        # processing that generates the \prad_tot signal in the tree (getbolo.pro
        # and powers.pro).  I converted them into Matlab routines, and modified the
        # analysis so that the smoothing is causal, and uses a shorter window.
        smoothing_window = 0.010  # [s]
        bol_prm = self.conn.get(r"\bol_prm").data()
        lower_names = [f"\bol_u{i+1:02d}_v"for i in range(24)]
        upper_names = [f"\bol_l{i+1:02d}_v" for i in range(24)]
        bol_chan_names = lower_names + upper_names
        bol_channels = []
        bol_top = []
        bol_time = []
        for i in range(48):
            bol_channels.append()
        bol_top = []
        for i in range(48):
            bol_time = self.conn.get(
                f"dim_of(\top.raw:bolo{i+1:02d})").data()/1.e3
        a_struct = get_bolo(self._shot_id, bol_channels,
                            bol_prm, bol_top, bol_time)
        ier = 0

    def get_ohmic_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        # Get edge loop voltage and smooth it a bit with a median filter
        try:
            v_loop = self.conn.get(
                f"ptdata('vloopb', {self._shot_id})").data()  # [V]
            t_v_loop = self.conn.get(
                f"dim_of(ptdata('vloopb', {self._shot_id})").data()/1.e3  # [ms]->[s]
            v_loop = scipy.signal.medfilt(v_loop, 11)
            v_loop = interp1(t_v_loop, v_loop, self._times, 'linear')
        except mdsExceptions.TreeFOPENR as e:
            v_loop = np.full(len(self._times), np.nan)
            t_v_loop = v_loop.copy()
       # Get plasma current
        try:
            ip = self.conn.get(f"ptdata('ip', {self._shot_id})").data()  # [A]
            t_ip = self.conn.get(
                f"dim_of(ptdata('ip', {self._shot_id})").data()/1.e3  # [ms] -> [s]
            # We choose a 20-point width for gsastd. This means a 10ms window for ip smoothing
            diptdt_smoothed = gsastd(t_ip, ip, 1, 20, 3, 1, 0)
            self.conn.openTree(self.efit_tree_name, self.shot_id)
            li = self.conn.get(r"\efit_a_eqdsk:li").data()
            t_li = self.conn.get(
                r"dim_of(\efit_a_eqdsk:li)").data()/1.e3  # [ms] -> [s]
            chisq = self.conn.get(r"\efit_a_eqdsk:chisq").data()
            # Filter out invalid indices of efit reconstruction
            invalid_indices = np.where
        except mdsExceptions.TreeFOPENR as e:
            p_ohm = np.full(len(self._times), np.nan)
            return p_ohm, v_loop
        # [m] For simplicity, use fixed r_0 = 1.67 for DIII-D major radius
        r_0 = 1.67
        inductance = 4.*np.pi*r_0 * li/2  # [H]
        inductance = interp1(t_li, inductance, self._times, 'linear')
        ip = interp1(t_ip, ip, self._times, 'linear')
        dipdt_smoothed = interp1(t_ip, dipdt_smoothed, self._times, 'linear')
        v_inductive = inductance * dipdt_smoothed  # [V]
        v_resistive = v_loop - v_inductive  # [V]
        p_ohm = ip * v_resistive  # [W]
        return p_ohm, v_loop

    # TODO: Cover all matlab density files

    def get_density_parameters(self):
        ne = np.full(len(self._times), np.nan)
        g_f = ne.copy()
        dne_dt = ne.copy()
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        try:
            t_ne = self.conn.get("dim_of(\density)").data()/1.e3  # [ms] -> [s]
            ne = self.conn.get("\density").data()*1.e6  # [cm^3] -> [m^3]
            dne_dt = np.gradient(ne, t_ne)
            ne = interp1(t_ne, ne, self._times, 'linear')
            dne_dt = interp1(t_ne, dne_dt, self._times, 'linear')
            t_ip = self.conn.get(
                f"dim_of(ptdata('ip', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip = self.conn.get(f"ptdata('ip', {self._shot_id})").data()  # [A]
            ipsign = np.sign(np.sum(ip))
            ip = interp1(t_ip, ip*ipsign, self._times, 'linear')
            a_minor = self.conn.get("\efit_a_eqdsk:aminor").data()  # [m]
            t_a = self.conn.get(
                "\efit_a_eqdsk:atime").data()/1.e3  # [ms] -> [s]
            a_minor = interp1(t_a, a_minor, self._times, 'linear')
            n_g = ip/1.e6 / (np.pi*a_minor**2)  # [MA/m^2]
            g_f = ne/1.e20 / n_g  # TODO: Fill in units
        except mdsExceptions.TreeFOPENR as e:
            # TODO: Better exception message
            # TODO: Confirm that there is a separate exception if ptdata name doesn't exist
            print("Failed to get some parameter")
        return ne, g_f, dne_dt

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
            a_minor_rt = self.conn.get("\efit_a_eqdsk:aminor").data()  # [m]
            t_a_rt = self.conn.get(
                "\efit_a_eqdsk:atime").data()/1.e3  # [ms] -> [s]
            a_minor_rt = interp1(t_a_rt, a_minor_rt, self._times, 'linear')
            n_g_rt = ip/1.e6 / (np.pi*a_minor_rt**2)  # [MA/m^2]
            g_f_rt = ne_rt/1.e20 / n_g_rt  # TODO: Fill in units
        except mdsExceptions.TreeFOPENR as e:
            # TODO: Better exception message
            print("Failed to get some parameter")
        return ne_rt, g_f_rt, dne_dt_rt

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
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get measured plasma current parameters")
        # Get programmed plasma current parameters
        try:
            t_ip_prog = self.conn.get(
                f"dim_of(ptdata('iptipp', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_prog = self.conn.get(
                f"ptdata('iptipp', {self._shot_id})").data()  # [A]
            polarity = np.unique(self.conn.get(
                f"ptdata('iptdirect', {self._shot_id})").data())
            if len(polarity) > 1:
                print("Polarity of Ip target is not constant")
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, self._times, 'linear')
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, self._times, 'linear')
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get programmed plasma current parameters")
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
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get ipimode signal")
            ipimode = np.full(len(self._times), np.nan)
        feedback_on_indices = np.where(ipimode == 0 | ipimode == 3)
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
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get epsoff signal")
            power_supply_railed = np.full(len(self._times), np.nan)
        return ip, ip_prog, ip_error, dip_dt, dipprog_dt, power_supply_railed

    def get_rt_ip_parameters(self):
        self.conn.openTree('d3d', self._shot_id)
        ip_rt = np.full(len(self._times), np.nan)
        ip_prog_rt = np.full(len(self._times), np.nan)
        ip_error_rt = np.full(len(self._times), np.nan)
        dip_dt_rt = np.full(len(self._times), np.nan)
        dipprog_dt_rt = np.full(len(self._times), np.nan)
        # Get measured plasma current parameters
        try:  # TODO: Ask about using ipspr15V
            t_ip_rt = self.conn.get(
                f"dim_of(ptdata('ipsip', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_rt = self.conn.get(
                f"ptdata('ipsip', {self._shot_id})").data()*1.e6  # [MA] -> [A]
            dip_dt_rt = np.gradient(ip_rt, t_ip_rt)
            ip_rt = interp1(t_ip_rt, ip_rt, self._times, 'linear')
            dip_dt_rt = interp1(t_ip_rt, dip_dt_rt, self._times, 'linear')
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get measured plasma current parameters")
        # Get programmed plasma current parameters
        try:
            t_ip_prog_rt = self.conn.get(
                f"dim_of(ptdata('ipsiptargt', {self._shot_id}))").data()/1.e3  # [ms] -> [s]
            ip_prog_rt = self.conn.get(
                f"ptdata('ipsiptargt', {self._shot_id})").data()*1.e6*.5  # [MA] -> [A]
            polarity = np.unique(self.conn.get(
                f"ptdata('iptdirect', {self._shot_id})").data())
            if len(polarity) > 1:
                print("Polarity of Ip target is not constant")
                polarity = polarity[0]
            ip_prog_rt = ip_prog_rt * polarity
            dipprog_dt_rt = np.gradient(ip_prog_rt, t_ip_prog_rt)
            ip_prog_rt = interp1(t_ip_prog_rt, ip_prog_rt,
                                 self._times, 'linear')
            dipprog_dt_rt = interp1(
                t_ip_prog_rt, dipprog_dt_rt, self._times, 'linear')
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get programmed plasma current parameters")
        try:
            t_ip_error_rt = self.conn.get(
                f"dim_of(ptdata('ipeecoil', {self._shot_id}))").data()/1.e3  # [ms] to [s]
            ip_error_rt = self.conn.get(
                f"ptdata('ipeecoil', {self._shot_id})").data()*1.e6*.5  # [MA] -> [A]
            ip_error_rt = interp1(
                t_ip_error_rt, ip_error_rt, self._times, 'linear')
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get ipeecoil signal")
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
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get ipimode signal")
            ipimode = np.full(len(self._times), np.nan)
        feedback_off_indices = np.where(ipimode != 0 & ipimode == 3)
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
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get epsoff signal")
            power_supply_railed = np.full(len(self._times), np.nan)
        return ip_rt, ip_prog_rt, ip_error_rt, dip_dt_rt, dipprog_dt_rt, power_supply_railed

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
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get vpszp signal")
            return z_cur, z_cur_norm, z_prog, z_error, z_error_norm
        self.conn.openTree(self.efit_tree_name, self._shot_id)
        try:
            t_a = self.conn.get(
                r"\efit_a_eqdsk:atime").data()/1.e3  # [ms] -> [s]
            a_minor = self.conn.get(r"\efit_a_eqdsk:aminor").data()  # [m]
            chisq = self.conn.get(r"\efit_a_eqdsk:chisq").data()
            invalid_indices = np.where(chisq > 50)
            a_minor[invalid_indices] = np.nan
            a_minor = interp1(t_a, a_minor, self._times, 'linear')
            z_cur_norm = z_cur/a_minor
        except mdsExceptions.TreeFOPENR as e:
            print("Failed to get efit parameters")
            z_cur_norm = z_cur / self.nominal_flattop_radius
        return z_cur, z_cur_norm, z_prog, z_error, z_error_norm

    # TODO: Complete n1 bradial method
    def get_n1_bradial(self):
        pass


if __name__ == '__main__':
    pass
    # shot = D3DShot('', D3D_DISRUPTED_SHOT)
    # print(shot.data.head())
