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
    pass 
    
class CmodShot(Shot):
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
        self.data['ip'],self.data['dip'],self.data['dip_smoothed'],self.data['ip_prog'],self.data['dipprog_dt'],self.data['ip_error'] = self._calc_Ip_parameters()
        self.data['p_oh'], self.data['v_loop'] = self._calc_p_ohm_v_loop()
        self.data['p_rad'], self.data['dprad_dt'], self.data['p_lh'], self.data['p_icrf'], self.data['p_input'], self.data['radiated_fraction'] = self._calc_power()
        self.data['kappa_area'] = self._calc_kappa_area()
        self.data['v_0'] = self._calc_rotation_velocity()
        self.data['v_0_uncalibrated'] = [np.nan]*len(self.data) #TODO: Populate when shot is missing from calibrated.txt 
        self.data['v_mid'] = [np.nan]*len(self.data) #TODO: Ask about removing from database
        self.data['sxr'] = self._calc_sxr_data()
        self.data['beta_n'], self.data['beta_p'], self.data['dbetap_dt'], self.data['kappa'], self.data['upper_gap'], self.data['lower_gap'], self.data['li'], self.data['dli_dt'], self.data['q0'], self.data['qstar'], self.data['q95'], _, self.data['Wmhd'], self.data['dWmhd_dt'], self.data['ssep'], self.data['n_over_ncrit'] = self._calc_EFIT_parameters()   
        self.data['z_error'], self.data['z_prog'], self.data['z_cur'], self.data['v_z'], self.data['z_times_v_z'] = self._calc_Z_parameters()
        self.data['iefc'] = self._calc_efc_current()
        self.data['n_equal_1_mode'],self.data['n_equal_1_normalized'], _ = self._calc_n_equal_1_amplitude()
        # Check there are no floats(float32) since the SQL database only contains type doubles(float64)
        float_columns = list(self.data.select_dtypes(include=['float32']).columns)
        print(float_columns) #TODO: Change to an assert
        # self.data[float_columns] = self.data[float_columns].astype('float64')

    def get_active_wire_segments(self):
        pcs_tree = Tree('pcs',self._shot_id)
        segment_nodes = pcs_tree.getNodeWild("\\top.seg_*")
        # Collect active segments and their information
        active_segments = []
        for node in segment_nodes: 
            if node.isOn():
                active_segments.append([node,node.getNode(":start_time").getData().data()])
        active_segments.sort(key=lambda n: n[1])
        end_times = np.roll(np.asarray([n[1] for n in active_segments]),-1)
        end_times[-1] = 12.383
        for i in range(len(active_segments)):
            active_segments[i].append(end_times[i])
        return active_segments

    @staticmethod
    def calc_time_until_disrupt():
        pass 

    def _calc_time_until_disrupt(self):
        if self._metadata['disrupted']:
            if self._metadata['disrupted'] == -1:
                raise ValueError("Shot disrupted but no t_disrupt found.")   
            time_until_disrupt = self._metadata['disrupted'] - self.data['times'] 
        return time_until_disrupt

    @staticmethod 
    def calc_IP_parameters(times, ip, magtime,ip_prog, pcstime):
        dip = np.gradient(ip, magtime)
        dip_smoothed = smooth(dip, 11)
        dipprog_dt = np.gradient(ip_prog,pcstime)
        ip_prog = interp1(pcstime, ip_prog,times, fill_value=ip_prog[-1])
        dipprog_dt = interp1(pcstime, dipprog_dt, times)
        ip = interp1(magtime, ip, times)
        dip = interp1(magtime, dip, times)
        dip_smoothed = interp1(magtime, dip_smoothed, times)
        ip_error = ip-ip_prog
        return ip, dip, dip_smoothed, ip_prog,dipprog_dt, ip_error

    def _calc_Ip_parameters(self):
        magnetics_tree = Tree('magnetics',self._shot_id) # Automatically generated
        active_segments = self.get_active_wire_segments()
        # Default PCS timebase is 1 KHZ
        pcstime = np.array(np.arange(-4,12.383,.001))       
        ip_prog = np.empty(pcstime.shape)
        ip_prog.fill(np.nan)
        # For each activate segment:
        # 1.) Find the wire for IP control and check if it has non-zero PID gains
        # 2.) IF it does, interpolate IP programming onto the PCS timebase 
        # 3.) Clip to the start and stop times of PCS timebase
        for segment,start,end in active_segments:
            # Ip wire can be one of 16 but is normally no. 16
            for wire_index in range(16,0,-1):
                wire_node = segment.getNode(f":P_{wire_index :02d}:name")
                if wire_node.getData().data() == 'IP':
                    try:
                        pid_gains = wire_node.getNode(":pid_gains").getData().data()
                        if np.any(pid_gains):
                            sig_node = segment.getNode(f":P_{wire_index :02d}")
                            signal_record = sig_node.getData()
                            sigtime = signal_record.dim_of(0)
                            signal = signal_record.data()
                            ip_prog = interp1(sigtime, signal, pcstime)
                            ip_prog = ip_prog[np.where(pcstime >= start &  pcstime <= end)]
                    except mdsExceptions.MdsException as e:
                        pass #TODO: Change
                else:
                    continue 
                break
        ip = magnetics_tree.getNode(r"\ip").getData().data().astype('float64',copy=False)
        magtime = magnetics_tree.getNode(r"\ip").getData().dim_of(0)
        return Shot.calc_IP_parameters(self._times, ip, magtime, ip_prog, pcstime)
    
    @staticmethod 
    def calc_Z_parameters(times, z_prog, pcstime, z_error_without_ip, ip, dpcstime):
        z_error = z_error_without_ip/ip # [m]
        z_prog_dpcs = interp1(pcstime, z_prog, dpcstime)
        z_cur = z_prog_dpcs + z_error #[m]
        v_z = np.gradient(z_cur, dpcstime) # m/s
        z_times_v_z = z_cur * v_z # m^2/s
        z_prog = interp1(pcstime, z_prog, times, 'linear', False, z_prog[-1])
        z_error = -interp1(dpcstime, z_error, times, 'linear', False,z_error[-1])
        z_cur = -interp1(dpcstime, z_cur, times, 'linear',False, z_cur[-1])
        v_z = interp1(dpcstime, v_z, times, 'linear',False, v_z[-1])
        z_times_v_z = interp1(dpcstime, z_times_v_z, times, 'linear', False, z_times_v_z[-1])
        return z_error, z_prog, z_cur, v_z, z_times_v_z

    def _calc_Z_parameters(self):
        pcstime = np.array(np.arange(-4,12.383,.001))
        z_prog = np.empty(pcstime.shape)
        z_prog.fill(np.nan)
        z_prog_temp = z_prog.copy()
        z_wire_index = -1
        active_wire_segments = self.get_active_wire_segments()
        for segment,start,end in active_wire_segments:
            for wire_index in range(1, 17):
                wire_node = segment.getNode(f":P_{wire_index :02d}:name")
                if wire_node.getData().data() == "ZCUR":
                    try:
                        pid_gains = segment.getNode(f":P_{wire_index :02d}:pid_gains").getData().data()
                        if np.any(pid_gains):
                            sig_node = segment.getNode(f":P_{wire_index :02d}")
                            signal_record = sig_node.getData()
                            sigtime = signal_record.dim_of(0)
                            signal = signal_record.data()
                            z_prog_temp = interp1(sigtime, signal, pcstime,'linear',False, fill_value = signal[-1])
                            z_wire_index = wire_index
                            segment_indices = [np.where((pcstime >= start) &  (pcstime <= end))]
                            z_prog[segment_indices] = z_prog_temp[segment_indices]
                            break
                    except mdsExceptions.MdsException as e:
                        print(e)
                        continue #TODO: Consider raising appropriate error
                else:
                    continue 
                break
        if z_wire_index == -1:
            raise ValueError("No ZCUR wire was found") #TODO: Make appropriate error
        # Read in A_OUT, which is a 16xN matrix of the errors for *all* 16 wires for
        # *all* of the segments. Note that DPCS time is usually taken at 10kHz. 
        hybrid_tree = Tree('hybrid', self._shot_id)
        wire_errors_record = hybrid_tree.getNode(r'\top.hardware.dpcs.signals:a_out').getData()
        wire_errors, dpcstime = wire_errors_record.data(), np.array(wire_errors_record.dim_of(1)) # s
        # The value of Z_error we read is not in the units we want. It must be *divided* by a factor AND *divided* by the plasma current.
        z_error_without_factor_and_ip = wire_errors[:,z_wire_index]
        z_error_without_ip = np.empty(z_error_without_factor_and_ip.shape)
        z_error_without_ip.fill(np.nan)
        # Also, it turns out that different segments have different factors. So we
        # search through the active segments (determined above), find the factors,
        # and *divide* by the factor only for the times in the active segment (as
        # determined from start_times and stop_times.
        for i in range(len(active_wire_segments)):
            segment, start,end = active_wire_segments[i]
            z_factor = hybrid_tree.getNode(fr'\dpcs::top.seg_{i+1:02d}:p_{z_wire_index:02d}:predictor:factor').getData().data()
            z_error_without_ip[np.where((dpcstime >= start) & (dpcstime <= end))] /= z_factor # [A*m] 
        #Next we grab ip, which comes from a_in:input_056. This also requires 
        #*multiplication* by a factor.
        #NOTE that I can't get the following ip_without_factor to work for shots
        #before 2015. 
        #TODO: Try to fix this
        if self._shot_id > 1150101000:
            ip_without_factor = hybrid_tree.getNode(r'\hybrid::top.hardware.dpcs.signals.a_in:input_056').getData().data()
            ip_factor = hybrid_tree.getNode(r'\hybrid::top.dpcs_config.inputs:input_056:p_to_v_expr').getData().data()
            ip = ip_without_factor*ip_factor #[A]
        else:
            magnetics_tree = Tree('magnetics',self._shot_id)
            ip_record = magnetics_tree.getNode('\ip').getData()
            ip = ip_record.data()
            ip_time = ip_record.dim_of(0)
            ip = interp1(ip_time,ip,dpcstime)
        return self.calc_Z_parameters(self._times,z_prog,pcstime, z_error_without_ip, ip, dpcstime)

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

    @staticmethod
    def calc_rotation_velocity(intensity, time, vel, hirextime):
        """ 
        Uses spectroscopy graphs of ionized(to hydrogen and helium levels) Argon to calculate velocity. Because of the heat profile of the plasma, suitable measurements are only found near the center
        """
        v_0 = np.empty(len(self._times))
        # Check that the argon intensity pulse has a minimum count and duration threshold
        valid_indices = np.where(intensity>1000 & intensity<10000)
        # Matlab code just multiplies by time delta but that doesn't work in the case where we have different time deltas
        # Instead we sum the time deltas for all valid indices to check the total duration
        if np.sum(time[valid_indices+1] - time[valid_indices]) >= .2:
            v_0 = interp1(hirextime, vel, self._times)
            v_0[np.where(abs(v_0) > 200)] = np.nan #TODO: Determine better threshold
            v_0 *= 1000.0
        return v_0

    # TODO: Calculate v_mid
    def _calc_rotation_velocity(self):
        calibrated = pd.read_csv("lock_mode_calib_shots.txt") #TODO: productionize
        # Check to see if shot was done on a day where there was a locked
        # mode HIREX calibration by cross checking with list of calibrated 
        # runs. If not calibrated, return NaN outputs.
        if self._shot_id not in calibrated:
            v_0 = np.empty(len(self._times))
            v_0.fill(np.nan)
            return v_0,v_mid
        try:
            spec_tree = Tree('spectroscopy', self._shot_id)
            intensity_record = spec_tree.getNode('.hirex_sr.analysis.a:int').getData()
            intensity = intensity_record.data().astype('float64',copy=False)
            time = intensity_record.dim_of(0)
            vel_record = spec_tree.getNode('.hirex_sr.analysis.a:vel').getData()
            vel = vel_record.data().astype('float64',copy=False)
            hirextime = vel_record.dim_of(0)
        except mdsExceptions.TreeFOPENR as e:
            print("WARNING: failed to open necessary tress for rotational velocity calculations.")
            v_0 = np.empty(len(self._times))
            v_0.fill(np.nan)
            return v_0
        return Shot.calc_rotation_velocity(intensity, time, vel, hirextime)

    # TODO: Split into static and instance method
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
        the n=1 signal. 

        N=1 toroidal assymmetry in the magnetic fields
        """
        n_equal_1_amplitude = np.empty(len(self._times))
        n_equal_1_amplitude.fill(np.nan)
        n_equal_1_normalized = n_equal_1_amplitude.copy()
        n_equal_1_phase = n_equal1_amplitude.copy()
        bp13_names = ['BP13_BC', 'BP13_DE', 'BP13_GH', 'BP13_JK'] # These sensors are placed toroidally around the machine. Letters refer to the 2 ports the sensors were placed between.
        bp13_signals = np.empty((len(self._times), len(bp13_names)))
        mag_tree = Tree('magnetics',self._shot_id)
        path = r"\mag_bp_coils"
        bp_node_names = mag_tree.getNode(path + "nodename").getData().data()
        phi = mag_tree.GetNode(path + "btor_pickup").getData().data()
        bp13_indices = np.intersect1d(bp_node_names,bp13_names)
        bp13_phi = phi[bp13_indices] + 360 # INFO
        bp13_btor_pickup_coeffs = btor_pickup_coeffs[bp13_indices]
        btor_record = mag_tree.getNode(path + r"\btor").getData()
        btor = btor_record.data()
        t_mag = btor_record.dim_of(0)
        baseline_indices = np.where(t_mag <= -1.8) # Toroidal power supply takes time to turn on, from ~ -1.8 and should be on by t=-1. So pick the time before that to calculate baseline
        btor = btor - np.mean(btor[baseline_indices])
        path = r"\mag_bp_coils.signals"
        # For each sensor:
        # 1. Subtract baseline offset
        # 2. Subtract btor pickup
        # 3. Interpolate bp onto shot timebase
        for i in range(len(bp13_names)):
            signal = mag_tree.getNode(path + bp13_names[i]).getData().data()
            if len(signal) == 1:
                print("WARNING: Can't fit with signal. Returning nans")
                return n_equal_1_amplitude, n_equal_1_normalized, n_equal_1_phase
            baseline = np.mean(signal[baseline_indices])
            signal = signal - baseline 
            signal = signal - bp13_btor_pickup_coeffs[i]*btor 
            bp13_signals[:,i] = interp1(t_mag, signal, self._times)
        polarity = np.sign() #TODO: Examine edge case behavior of sign 
        btor_magnitude = btor*polarity 
        btor_magnitude = interp1(t_mag, btor_magnitude, self._times)
        # Create the 'design' matrix ('A') for the linear system of equations:
        # Bp(phi) = A1 + A2*sin(phi) + A3*cos(phi)
        ncoeffs = 3 
        A = np.empty((len(bp13_names),ncoeffs))
        A[:,0] = np.ones
        A[:, 1] = np.sin(bp13_phi*np.pi/180.0)
        A[:,2] = np.cos(bp13_phi*np.pi/180.0)
        coeffs = np.linalg.inv(A) @ bp13_signals 
        # The n=1 amplitude at each time is sqrt(A2^2 + A3^2)
        # The n=1 phase at each time is arctan(-A2/A3), using complex number
        # phasor formalism, exp(i(phi - delta))
        n_equal_1_amplitude = np.sqrt(coeffs[1,:]**2 + coeffs[2,:]**2)
        n_equal_1_phase = np.arctan2(-coeffs[1,:], coeffs[2,:]) #TODO: Confirm arctan2 = atan2
        n_equal_1_normalized = n_equal_1_ampltiude / btor_magnitude

        # INFO: Debugging purpose block of code at end of matlab file
        # INFO: n_equal_1_amplitude vs n_equal_1_mode
        return n_equal_1_amplitude, n_equal_1_normalized, n_equal_1_phase




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
        
    @staticmethod 
    def calc_efc_current(times, iefc, t_iefc):
        return interp1(t_iefc, iefc, times, 'linear')

    def _calc_efc_current(self):
        eng_tree = Tree('engineering',self._shot_id)
        iefc_record = eng_tree.getNode(r'\efc:_u_bus_r_cur').getData()
        iefc, t_iefc = iefc_record.data(), iefc_record.dim_of(0)
        return calc_efc_current(self._times,iefc, t_iefc)

    #TODO: Split 
    @staticmethod 
    def calc_Ts_data(times, ts_data, ts_time, ts_z):
        zarray = np.array(np.arange(-.5,.5,.001))
        valid_times = np.where(ts_time > 0)
        y = ts_data[valid_times,:]
        ok_indices = np.where(y != 0)
        y = y[ok_indices]
        z = ts_z[ok_indices]
        _, _, sigma = gaussian_fit(z,y)
        te_hwm = sigma*1.1774 #50%
        te_hwm = interp1(ts_time, te_hwm, self._times)
        return te_hwm 

    def _calc_Ts_data(self):
        # TODO: Guassian vs parabolic fit for te profile
        te_hwm = np.empty((len(self._times)))
        electron_tree = Tree("electrons",self._shot_id)

        # Read in Thomson core temperature data, which is a 2-D array, with the
        # dependent dimensions being time and z (vertical coordinate)
        node_path = ".yag_new.results.profiles"
        try:
            ts_data = electron_tree.getNode(node_path + ":te_rz").getData().data()
            ts_time = electron_tree.getNode(node_path + ":te_rz").getData().dim_of(0)
            ts_z = electron_tree.getNode(node_path + ":z_sorted").getData().data()
        except mdsExceptions.MdsException as e:
            print(e) #TODO: Change 
            te_hwm.fill(np.nan)
            return te_hwm 
        return Shot.calc_Ts_data(self._times, ts_data, ts_time, ts_z)

    # TODO: Finish
    @staticmethod 
    def calc_peaking_factor():
        pass 

    def _calc_peaking_factor(self):
        pass        

    @staticmethod 
    def calc_sxr_data():
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
    efit_vars = {'beta_n':'\efit_a_eqdsk:betan', 'beta_p':'\efit_a_eqdsk:betap', 'kappa':'\efit_a_eqdsk:kappa', 'li':'\efit_a_eqdsk:li', 'upper_gap':'\efit_a_eqdsk:gaptop', 'lower_gap':'\efit_a_eqdsk:gapbot',  'q0':'\efit_a_eqdsk:q0', 'qstar': '\efit_a_eqdsk:qstar','q95':'\efit_a_eqdsk:q95',  'v_loop_efit':'\efit_a_eqdsk:vsurf', 'Whmd':'\efit_a_eqdsk:wmhd', 'chisq':'\efit_a_eqdsk:chisq'}
    efit_derivs = ['\efit_a_eqdsk:betap','\efit_a_eqdsk:li','Wmhd']
    def __init__(self, shot_id,efit_tree_name, data_columns = DEFAULT_SHOT_COLUMNS,data = None):
        self._shot_id = shot_id
        self._metadata = {
            'labels': {},
            'commit_hash': subprocess.check_output(["git", "describe","--always"]).strip(), 
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  #TODO: Fix
        }
        self.data = data
		self.efit_tree_name = efit_tree_name
        if data is None:
            self.data = pd.DataFrame()
            self._populate_shot_data()
	
    def _populate_shot_data(self):
        efit_data = self._calc_efit_parameters()
        self.data = pd.concat([self.data, efit_data], ignore_index=True)
	
    @staticmethod
    def calc_efit_parameters():
        pass
    def _calc_efit_parameters(self):
        efit_tree = Tree(self.efit_tree_name,self._shot_id)
        efit_data = {k:efit_tree.getNode(v) for k,v in self.efit_vars}
        efit_data['time']/=1000
        for key in efit_derivs:
            efit_data['d' + key + '_dt'] = np.gradient(efit_data[key],efit_data['time']
        return pd.DataFrame([efit_data]) 
        		



class EASTShot(Shot):
    pass
