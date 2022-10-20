import math
import subprocess 

import pandas as pd
import numpy as np 
from scipy.interpolate import interp1d
from scipy.signal.windows import boxcar

from MDSplus import *

from utils import interp1, smooth

DEFAULT_SHOT_COLUMNS = ['time','shot','time_until_disrupt','ip']

class Shot:
    """
    Base shot class to represent a single shot of a fusion experiment.

    Attributes
    ---
    metadata->dict: 
        labels: {} /* Dictionary where key is the label name(e.g. "disruption period")
                    and value is a list of starting and ending time tuples(e.g.[(10,100)]*/
        commit_hash: // Commit hash of code version that was used to generate shot data 
        timestep: // Time delta. Every array must conform to this time step
        duration: // Length in milliseconds of shot
        description: // Optional text summary to add to shot
    data->pandas.DataFrame:
    """

    # TODO: Populate metadata dict
    def __init__(self, mdsplus_name,shot_id,data_columns = DEFAULT_SHOT_COLUMNS,data = None):
        # self._tree = Tree(mdsplus_name,shot_id)
        self._shot_id = shot_id
        self._metadata = {
            'labels': {},
            'commit_hash': subprocess.check_output(["git", "describe","--always"]).strip(), 
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  #TODO: Fix
        }
        # Analysis tree has automatically generated Efit data. Efi18 is manually created Efit data
        try:
            self._analysis_tree = Tree('efit18', shot_id, mode="readonly")
            self._times = self._analysis_tree.getNode(r"\efit18::efit_aeqdsk:time").getData().data().astype('float64',copy=False)
        except mdsExceptions.TreeFOPENR as e:
            try:
                self._analysis_tree = Tree('analysis', shot_id, mode="readonly")
                self._times = self._analysis_tree.getNode(r"\analysis::efit_aeqdsk:time").getData().data().astype('float64',copy=False)
            except mdsExceptions.TreeFOPENR as f: 
                print("WARNING: No EFIT data found")
        self.data = data
        if data is None:
            self.data = pd.DataFrame()
            self._populate_shot_data()
    
    def _populate_shot_data(self):
        self.data['times'] = self._times
        self.data['commit_hash'] = self._metadata['commit_hash'] #TODO: convert from bytes
        self.data['time_until_disrupt'] = self._calc_time_until_disrupt() 
        self.data['ip'],self.data['dip'],self.data['dip_smoothed'] = self._calc_Ip_parameters()
        self.data['p_ohm'], self.data['v_loop'] = self._calc_p_ohm_v_loop()
        self.data['p_rad'], self.data['dprad'], self.data['p_lh'], self.data['p_icrf'], self.data['p_input'], self.data['rad_fraction'] = self._calc_power()
        self.data['kappa_area'] = self._calc_kappa_area()
        self.data['v_0'],self.data['v_mid'] = self._calc_rotation_velocity()
        self.data['sxr'] = self._calc_sxr_data()
        self.data['beta_N'], self.data['beta_p'], self.data['beta_p_dot'], self.data['kappa'], self.data['upper_gap'], self.data['lower_gap'], self.data['li'], self.data['li_dot'], self.data['q0'], self.data['qstar'], self.data['q95'], self.data['V_loop_efit'], self.data['Wmhd'], self.data['dWmhd_dt'], self.data['ssep'], self.data['n_over_ncrit'] = self._calc_EFIT_parameters()   
        # Check there are no floats(float32) since the SQL database only contains type doubles(float64)
        float_columns = list(self.data.select_dtypes(include=['float32']).columns)
        print(float_columns)
        # self.data[float_columns] = self.data[float_columns].astype('float64')

    def _calc_time_until_disrupt(self):
        if self._metadata['disrupted']:
            if self._metadata['disrupted'] == -1:
                raise ValueError("Shot disrupted but no t_disrupt found.")   
            time_until_disrupt = self._metadata['disrupted'] - self.data['times'] 
        return time_until_disrupt

    # TODO: Fix plasma current calculation. Add programmed current
    @staticmethod 
    def calc_IP_parameters(times, ip, magtime):
        dip = np.gradient(ip, magtime)
        dip_smoothed = smooth(dip, 11)
        ip = interp1(magtime, ip, times)
        dip = interp1(magtime, dip, times)
        dip_smoothed = interp1(magtime, dip_smoothed, times)
        return ip, dip, dip_smoothed

    def _calc_Ip_parameters(self):
        magnetics_tree = Tree('magnetics',self._shot_id) # Automatically generated
        ip = magnetics_tree.getNode(r"\ip").getData().data().astype('float64',copy=False)
        magtime = magnetics_tree.getNode(r"\ip").getData().dim_of(0)
        return Shot.calc_IP_parameters(self._times, ip, magtime)

    def _calc_Z_parameters(self):
        raise NotImplementedError("Skipping over this function for now because not core to demonstrating library capabilities")

    @staticmethod
    def calc_p_oh_v_loop(times, v_loop,v_loop_time,li,efittime, dip_smoothed, ip):
        inductance = 4.0*np.pi*1.0e-7 * 0.68 * li/2.0 # For simplicity, we use R0 = 0.68 m, but we could use \efit_aeqdsk:rmagx
        v_loop = interp1(v_loop_time, v_loop, times)
        inductance = interp1(efittime, inductance, times)
        v_inductive = inductance * dip_smoothed
        v_resistive = v_loop - v_inductive
        p_ohm = ip * v_resistive
        return p_ohm, v_loop       

    def _calc_p_ohm_v_loop(self):
        v_loop_record = self._analysis_tree.getNode(r"\top.mflux:v0").getData()
        v_loop = v_loop_record.data().astype('float64',copy=False)
        v_loop_time = v_loop_record.dim_of(0)
        if len(v_loop_time) <= 1:
            return None, None
        li_record = self._analysis_tree.getNode(r"\efit_aeqdsk:li").getData()
        li = li_record.data().astype('float64',copy=False)
        efittime = li_record.dim_of(0)
        return Shot.calc_p_oh_v_loop(self._times,v_loop, v_loop_time, li, efittime, self.data['dip_smoothed'],self.data['ip'])

    @staticmethod 
    def calc_power(times, p_lh, t_lh,p_icrf, t_icrf, p_rad, t_rad,p_ohm):
        p_lh = interp1(t_lh, p_lh * 1.0e3, times, bounds_error = False) if p_lh is not None else np.zeros(len(times))
        p_icrf = interp1(t_icrf, p_icrf * 1.0e6, times,bounds_error = False) if p_icrf is not None else np.zeros(len(times))
        if len(t_rad) == 1 or p_rad is None:
            p_rad = np.array([np.nan]*len(times)) #TODO: Fix 
            dprad = p_rad.copy()
        else: 
            dprad = np.gradient(p_rad,t_rad)
            p_rad = interp1(t_rad, p_rad, times)
            dprad = interp1(t_rad, dprad, times)   
        p_input = p_ohm+ p_lh + p_icrf 
        rad_fraction = p_rad/p_input
        rad_fraction[rad_fraction == np.inf] = np.nan 
        return p_rad, dprad, p_lh, p_icrf, p_input, rad_fraction

    def _calc_power(self):
        """
        NOTE: the timebase for the LH power signal does not extend over the full
            time span of the discharge.  Therefore, when interpolating the LH power
            signal onto the "timebase" array, the LH signal has to be extrapolated
            with zero values.  This is an option in the 'interp1' routine.  If the
            extrapolation is not done, then the 'interp1' routine will assign NaN
            (Not-a-Number) values for times outside the LH timebase, and the NaN's
            will propagate into p_input and rad_fraction, which is not desirable.
        """
        values = [None]*6
        trees = ['LH','RF','spectroscopy']
        nodes = [r'\LH::TOP.RESULTS:NETPOW',r"\rf::rf_power_net",r"\twopi_diode"]
        for i in range(3):
            try:
                tree = Tree(trees[i],self._shot_id)
                record = tree.getNode(nodes[i])
                values[2*i] = record.data().astype('float64',copy=False)
                values[2*i +1] = record.dim_of(0)
            except mdsExceptions.TreeFOPENR as e:
                continue
        return Shot.calc_power(self._times, *values, self.data['p_ohm'])

    def _calc_EFIT_parameters(self):
        efittime = self._analysis_tree.getNode('\efit_aeqdsk:time')
        beta_N = self._analysis_tree.getNode('\efit_aeqdsk:betan').getData().data().astype('float64',copy=False) #
        beta_p = self._analysis_tree.getNode('\efit_aeqdsk:betap').getData().data().astype('float64',copy=False) #
        kappa = self._analysis_tree.getNode('\efit_aeqdsk:eout').getData().data().astype('float64',copy=False) #
        li = self._analysis_tree.getNode('\efit_aeqdsk:li').getData().data().astype('float64',copy=False) #
        upper_gap = self._analysis_tree.getNode('\efit_aeqdsk:otop').getData().data().astype('float64',copy=False) #100 % meters
        lower_gap = self._analysis_tree.getNode('\efit_aeqdsk:obott').getData().data().astype('float64',copy=False) #100 % meters
        q0 = self._analysis_tree.getNode('\efit_aeqdsk:q0').getData().data().astype('float64',copy=False) #
        qstar = self._analysis_tree.getNode('\efit_aeqdsk:qstar').getData().data().astype('float64',copy=False) #
        q95 = self._analysis_tree.getNode('\efit_aeqdsk:q95').getData().data().astype('float64',copy=False) #
        V_loop_efit = self._analysis_tree.getNode('\efit_aeqdsk:vloopt').getData().data().astype('float64',copy=False) #
        Wmhd = self._analysis_tree.getNode('\efit_aeqdsk:wplasm').getData().data().astype('float64',copy=False) #
        ssep = self._analysis_tree.getNode('\efit_aeqdsk:ssep').getData().data().astype('float64',copy=False) #100 % meters
        n_over_ncrit =self._analysis_tree.getNode('\efit_aeqdsk:xnnc').getData().data().astype('float64',copy=False) #
        beta_p_dot = np.gradient(beta_p, efittime)
        li_dot = np.gradient(li, efittime)
        dWmhd_dt = np.gradient(Wmhd, efittime)
        beta_N = interp1(efittime, beta_N, self._times)
        beta_p = interp1(efittime, beta_p, self._times)
        kappa = interp1(efittime, kappa, self._times)
        li = interp1(efittime, li, self._times)
        upper_gap = interp1(efittime, upper_gap, self._times)
        lower_gap = interp1(efittime, lower_gap, self._times)
        q0= interp1(efittime, q0, self._times)
        qstar = interp1(efittime, qstar, self._times)
        q95 = interp1(efittime, q95, self._times)
        V_loop_efit = interp1(efittime, V_loop_efit, self._times)
        Wmhd = interp1(efittime, Wmhd, self._times)
        beta_p_dot = interp1(efittime, beta_p_dot, self._times)
        li_dot = interp1(efittime, li_dot, self._times)
        dWmhd_dt = interp1(efittime, dWmhd_dt, self._times)
        ssep = interp1(efittime, ssep, self._times)
        n_over_ncrit = interp1(efittime, n_over_ncrit, self._times)
        return beta_N, beta_p, beta_p_dot, kappa, upper_gap, lower_gap, li, li_dot, q0, qstar, q95, V_loop_efit, Wmhd, dWmhd_dt, ssep, n_over_ncrit

    @staticmethod 
    def calc_kappa_area(times, aminor, area,a_times):
        return interp1(a_times,area/(np.pi * aminor**2), times)

    def _calc_kappa_area(self):
        aminor = self._analysis_tree.getNode(r'\efit_a_eqdsk:aminor').getData().data().astype('float64',copy=False)
        area = self._analysis_tree.getNode(r'\efit_a_eqdsk:area').getData().data().astype('float64',copy=False)
        times = self._analysis_tree.getNode(r'\efit_aeqdsk:time').getData().data().astype('float64',copy=False)
        return Shot.calc_kappa_area(self._times, aminor, area, times)

    #TODO: Consider if necessary
    @staticmethod
    def calc_rotation_velocity():
        pass

    # TODO: Calculate v_mid
    def _calc_rotation_velocity(self):
        calibrated = pd.read_csv("lock_mode_calib_shots.txt") #TODO: productionize
        v_0 = np.empty(len(self._times))
        v_mid = np.empty(len(self._times))
        # Check to see if shot was done on a day where there was a locked
        # mode HIREX calibration by cross checking with list of calibrated 
        # runs. If not calibrated, return NaN outputs.
        if self._shot_id not in calibrated:
            v_0.fill(np.nan)
            v_mid.fill(np.nan)
            return v_0,v_mid
        try:
            spec_tree = Tree('spectroscopy', self._shot_id)
            intensity_record = spec_tree.getNode('.hirex_sr.analysis.a:int').getData()
            intensity = intensity_record.data().astype('float64',copy=False)
            time = intensity_record.dim_of(0)
            vel_record = spec_tree.getNode('.hirex_sr.analysis.a:vel').getData()
            vel = vel_record.data().astype('float64',copy=False)
            hirextime = vel_record.dim_of(0)
            # Check that the argon intensity pulse  a minimum count and duration threshold
            valid_indices = np.where(intensity>1000 & intensity<10000)
            # Matlab code just multiplies by time delta but that doesn't work in the case where we have different time deltas
            # Instead we sum the time deltas for all valid indices
            if np.sum(time[valid_indices+1] - time[valid_indices]) >= .2:
                set_v_0_interp = interp1d(hirextime, vel, kind='linear')
                v_0 = set_v_0_interp(self._times)
                v_0[np.where(abs(v_0) > 200)] = np.nan #TODO: Determine better threshold
                v_0 *= 1000.0
        except mdsExceptions.TreeFOPENR as e:
            print("WARNING: failed to open necessary tress for rotational velocity calculations.")
        return v_0,v_mid
        
    @staticmethod 
    def calc_n_equal_1_amplitude():
        pass 

    def _calc_n_equal_1_amplitude(self):
        """ Calculate n=1 amplitude and phase for the disruption warning database
        using the four BP13 Bp sensors near the midplane on the outboard vessel
        wall.  The calculation is done by using a least squares fit to an
        expansion in terms of n = 0 & 1 toroidal harmonics.  The BP13 sensors are
        part of the set used for plasma control and equilibrium reconstruction,
        and their signals have been analog integrated (units: tesla), so they
        don't have to be numerically integrated.  These four sensors were working
        well in 2014, 2015, and 2016.  I looked at our locked mode MGI run on
        1150605, and the different applied A-coil phasings do indeed show up on
        the n=1 signal. """
        pass

    @staticmethod 
    def calc_densities(times, n_e,t_n, ip, t_ip,a_minor,t_a):
        if len(n_e) == len(t_n):
            nan_arr = np.empty(len(times))
            nan_arr.fill(np.nan)
            return arr,arr,arr
        dn_dt = np.gradient(n_e,t_n)
        n_e = interp1(t_n,n_e,times)
        dn_dt = interp1(t_n,dn_dt,times)
        ip = -ip/1e6 # Convert from A to MA and take positive value
        ip = interp1(t_ip,ip,times)
        a_minor = interp1(t_a,aminor,times)
        n_G = ip/(np.pi*a_minor**2)*1e20 # Greenwald density in m ^-3
        g_f = abs(n_e/n_G)
        return n_e, dn_dt, g_f   

    def _calc_densities(self):
        e_tree = Tree('electroncs', self._shot_id)
        n_e_record = e_tree.getNode(r'.tci.results:nl_04/0.6').getData()
        n_e = record.data().astype('float64',copy=False)
        t_n = record.dim_of(0)
        mag_tree = Tree('magnetics',self._shot_id)
        ip_record = mag_tree.getNode(r'\ip').getData()
        ip = ip_record.data().astype('float64',copy=False)
        t_ip = ip_record.dim_of(0)
        a_tree = Tree('analysis',self._shot_id)
        a_minor_record = a_tree.getNode(r'.efit.results.a_eqdsk:aminor').getData()
        t_a = a_minor_record.dim_of(0)
        a_minor = a_minor_record.data().astype('float64',copy=False)
        return Shot.calc_densities(self._times, n_e, t_n, ip, t_ip, a_minor,t_a)
        

    def _calc_efc_current(self):
        pass

    # TODO: get more accurate description of soft x-ray data
    def _calc_sxr_data(self):
        """ """
        sxr = np.empty(len(self._times))
        try:
            tree = Tree('xtomo', self._shot_id)
            sxr_record = tree.getNode(r'\top.brightnesses.array_1:chord_16').getData()
            sxr = sxr_record.data().astype('float64',copy=False)    
            t_sxr = sxr_record.dim_of(0)
            sxr = interp1(t_sxr, sxr, self._times)
        except mdsExceptions.TreeFOPENR as e:
            print("WARNING: Failed to get SXR data")    
        return sxr

    def __getitem__(self, key):
        return self._metadata if key == 'metadata' else self.data[key]

class D3DShot(Shot):
    pass

class EASTShot(Shot):
    pass

class CMODShot(Shot):

    def _calc_TS_data_cmod(self):
        pass
    