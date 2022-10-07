import math
import subprocess 

import pandas as pd
import numpy as np 
from scipy.interpolate import interp1d
from scipy.signal.windows import boxcar

from MDSplus import *

from utils import interp1

DEFAULT_SHOT_COLUMNS = ['time','shot','time_until_disrupt','ip']

nan_arr = np.empty(len(self._times))
nan_arr.fill(np.nan)

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
            self._times = self._analysis_tree.getNode(r"\efit18::efit_aeqdsk:time").getData().data()
        except mdsExceptions.TreeFOPENR as e:
            try:
                self._analysis_tree = Tree('analysis', shot_id, mode="readonly")
                self._times = self._analysis_tree.getNode(r"\analysis::efit_aeqdsk:time").getData().data()
            except mdsExceptions.TreeFOPENR as f: 
                print("WARNING: No EFIT data found")
        self.data = data
        if data is None:
            self.data = pd.DataFrame()
            self._populate_shot_data()
    
    def _populate_shot_data(self):
        self.data['times'] = self._times
        self.data['time_until_disrupt'] = self._calc_time_until_disrupt() 
        self.data['ip'],self.data['dip'],self.data['dip_smoothed'] = self._calc_Ip_parameters()
        self.data['p_ohm'], self.data['v_loop'] = self._calc_p_ohm_v_loop()
        self.data['p_rad'], self.data['dprad'], self.data['p_lh'], self.data['p_icrf'], self.data['p_input'], self.data['rad_fraction'] = self._calc_power()
        self.data['kappa_area'] = self._calc_kappa_area()
        self.data['v_0'],self.data['v_mid'] = self._calc_rotation_velocity()
        self.data['sxr'] = self._calc_sxr_data()
        self.data['beta_N'], self.data['beta_p'], self.data['beta_p_dot'], self.data['kappa'], self.data['upper_gap'], self.data['lower_gap'], self.data['li'], self.data['li_dot'], self.data['q0'], self.data['qstar'], self.data['q95'], self.data['V_loop_efit'], self.data['Wmhd'], self.data['dWmhd_dt'], self.data['ssep'], self.data['n_over_ncrit'] = self._calc_EFIT_parameters()   


    def _calc_time_until_disrupt(self):
        if self._metadata['disrupted']:
            if self._metadata['disrupted'] == -1:
                raise ValueError("Shot disrupted but no t_disrupt found.")   
            time_until_disrupt = self._metadata['disrupted'] - self.data['times'] 
        return time_until_disrupt

    # TODO: Fix plasma current calculation. Add programmed current
    def _calc_Ip_parameters(self):
        magnetics_tree = Tree('magnetics',self._shot_id) # Automatically generated
        ip = magnetics_tree.getNode(r"\ip").getData().data()
        magtime = magnetics_tree.getNode(r"\ip").getData().dim_of(0)
        dip = np.gradient(ip, magtime)
        dip_smoothed = smooth(dip, 11)
        ip = interp1(magtime, ip, self._times)
        dip = interp1(matime, dip, self._times)
        dip_smoothed = interp1(magtime, dip_smoothed, self._times)
        return ip, dip, dip_smoothed

    def _calc_Z_parameters(self):
        pass

    def _calc_p_ohm_v_loop(self):
        v_loop_record = self._analysis_tree.getNode(r"\top.mflux:v0").getData()
        v_loop = v_loop_record.data()
        v_loop_time = v_loop_record.dim_of(0)
        if len(v_loop_time) <= 1:
            return nan_arr.copy(),nan_arr.copy()
        li_record = self._analysis_tree(r"\efit_aeqdsk:li").getData()
        li = li_record.data()
        efittime = li_record.dim_of(0)
        inductance = 4*pi*1.e-7 * 0.68 * li/2 # For simplicity, we use R0 = 0.68 m, but we could use \efit_aeqdsk:rmagx
        v_loop = interp1(v_loop_time, v_loop, self._times)
        inductance = interp1(efittime, inductance, self._times)
        v_inductive = inductance * self.data['dip_smoothed']
        v_resistive = v_loop - v_inductive
        p_ohm = self.data['ip'] * v_resistive
        return p_ohm, v_loop

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
        lh_tree = Tree('LH',self._shot_id)
        lh_record = lh_tree.getNode('\LH::TOP.RESULTS:NETPOW').getData()
        p_lh = 1e3 *lh_record.data()
        t_lh = lh_record.dim_of(0)
        p_lh = interp1(t_lh, p_lh, self._times, bounds_error = False)
        try:
            rf_tree = Tree('RF',self._shot_id)
            p_icrf_record = rf_tree.GetNode(r"(\rf::rf_power_net").getData()
            p_icrf = 1e6 * p_icrf_record.data()
            t_icrf = p_icrf_record.dim_of(0)
            p_icrf = interp1(t_icrf, p_icrf, self._times,bounds_error = False)
        except mdsExceptions.TreeFOPENR as e:
            print(e)
            p_icrf = np.zeros(len(self._times))
        spec_tree = Tree('spectroscopy', self._shot_id)
        p_rad_record = spec_tree.getNode(r"\twopi_diode").getData()
        p_rad = p_rad_record.data()
        t_rad = p_rad_record.dim_of(0)
        if len(t_rad) == 1:
            return nan_arr,nan_arr
        # TODO: 
        dprad = gradient(p_rad,t_rad)
        p_rad = interp1(t_rad, p_rad, self._times)
        dprad = interp1(t_rad, dprad, self._times)
        p_input = p_ohm + p_lh + p_icrf 
        rad_fraction = p_rad/p_input
        rad_fraction[rad_fraction == np.inf] = np.nan 
        return p_rad, dprad, p_lh, p_icrf, p_input, rad_fraction

    def _calc_EFIT_parameters(self):
        efittime = self._analysis_tree.getNode('\efit_aeqdsk:time')
        beta_N = self._analysis_tree.getNode('\efit_aeqdsk:betan').getData().data() #
        beta_p = self._analysis_tree.getNode('\efit_aeqdsk:betap').getData().data() #
        kappa = self._analysis_tree.getNode('\efit_aeqdsk:eout').getData().data() #
        li = self._analysis_tree.getNode('\efit_aeqdsk:li').getData().data() #
        upper_gap = self._analysis_tree.getNode('\efit_aeqdsk:otop').getData().data() #100 % meters
        lower_gap = self._analysis_tree.getNode('\efit_aeqdsk:obott').getData().data() #100 % meters
        q0 = self._analysis_tree.getNode('\efit_aeqdsk:q0').getData().data() #
        qstar = self._analysis_tree.getNode('\efit_aeqdsk:qstar').getData().data() #
        q95 = self._analysis_tree.getNode('\efit_aeqdsk:q95').getData().data() #
        V_loop_efit = self._analysis_tree.getNode('\efit_aeqdsk:vloopt').getData().data() #
        Wmhd = self._analysis_tree.getNode('\efit_aeqdsk:wplasm').getData().data() #
        ssep = self._analysis_tree.getNode('\efit_aeqdsk:ssep').getData().data() #100 % meters
        n_over_ncrit =self._analysis_tree.getNode('-\efit_aeqdsk:xnnc').getData().data() #
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

    def _calc_kappa_area(self):
        aminor = self._analysis_tree.getNode(r'\efit_a_eqdsk:aminor').getData().data()
        area = self._analysis_tree.getNode(r'\efit_a_eqdsk:area').getData().data()
        times = self._analysis_tree.getNode()
        kappa_area = area/(np.pi * aminor**2)
        set_kappa_interp = interp1d(self._times, kappa_area, kind='linear')
        kappa_area = set_kappa_interp(self._times) #Redundant for now I know, this is a placeholder
        kappa_area = interp1(sel)
        return kappa_area

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
            intensity = intensity_record.data()
            time = intensity_record.dim_of(0)
            vel_record = spec_tree.getNode('.hirex_sr.analysis.a:vel').getData()
            vel = vel_record.data()
            hirextime = vel_record.dim_of(0)
            # Check that the argon intensity pulse  a minimum count and duration threshold
            valid_indices = np.where(intensity>1000 & intensity<10000)
            # Matlab code just multiplies by time delta but that doesn't work in the case where we have different time deltas
            # Instead we sum the time deltas for all valid indices
            if np.sum(time[valid_indices+1] - time[valid_indices]) >= .2:
                set_v_0_interp = interp1d(hirextime, vel, kind='linear')
                v_0 = set_v_0_interp(self._times)
                v_0[np.where(abs(v_0) > 200)] = np.nan #TODO: Determine better threshold
                v_0 *= 1000
        except mdsExceptions.TreeFOPENR as e:
            print("WARNING: failed to open necessary tress for rotational velocity calculations.")
        return v_0,v_mid
        
    
    def _calc_n_equal_1_amplitude(self):
        pass

    def _calc_densities(self):
        e_tree = Tree('electroncs', self._shot_id)
        n_e_record = e_tree.getNode(r'.tci.results:nl_04/0.6').getData()
        n_e = record.data()
        t_n = record.dim_of(0)
        if len(n_e) == len(t_n):
            nan_arr = np.empty(len(self._times))
            nan_arr.fill(np.nan)
            return arr,arr,arr
        dn_dt = np.gradient(n_e,t_n)
        n_e = interp1(t_n,n_e,self._times)
        dn_dt = interp1(t_n,dn_dt,self._times)
        mag_tree = Tree('magnetics',self._shot_id)
        ip_record = mag_tree.getNode(r'\ip').getData()
        ip = -ip_record.data()/1e6 # Convert from A to MA and take positive value
        t_ip = ip_record.dim_of(0)
        a_tree = Tree('analysis',self._shot_id)
        a_minor_record = a_tree.getNode(r'.efit.results.a_eqdsk:aminor').getData()
        ip = interp1(t_ip,ip,self._times)
        a_minor = interp1(a_minor_record.dim_of(0),a_minor_record.data(),self._times)
        n_G = ip/(np.pi*a_minor**2)*1e20 # Greenwald density in m ^-3
        g_f = abs(n_e/n_G)
        return n_e, dn_dt, g_f
        

    def _calc_efc_current(self):
        pass

    # TODO: get more accurate description of soft x-ray data
    def _calc_sxr_data(self):
        """ """
        sxr = np.empty(len(self._times))
        try:
            tree = Tree('xtomo', self._shot_id)
            sxr_record = tree.getNode(r'\top.brightnesses.array_1:chord_16').getData()
            sxr = sxr_record.data()    
            sxr_time = sxr_record.dim_of(0)
            set_sxr_interp = interp1d(sxr_time,sxr,kind='linear')
            sxr = set_sxr_interp(self._times)
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
    