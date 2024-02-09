import traceback
import numpy as np

import pandas as pd
from disruption_py.databases.database import ShotDatabase
from disruption_py.mdsplus_integration.mds_connection import MDSConnection
from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.parameter_methods.cmod.basic_parameter_methods import BasicCmodRequests
from disruption_py.shots.shot_manager import ShotManager
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.utils import without_duplicates


class D3DShotManager(ShotManager):
    
    # Disruption Variables
    DT_BEFORE_DISRUPTION = 0.002
    DURATION_BEFORE_DISRUPTION = 0.10
    
    def shot_setup(
        self,
        shot_id : int,
        shot_settings : ShotSettings,
        **kwargs
    ) -> ShotProps:
        """
        Sets up the shot properties for cmod.
        """
        
        try:
            disruption_time=self.process_database.get_disruption_time(shot_id=shot_id)
        except Exception as e:
            disruption_time=None
            self.logger.error(f"Failed to retreive disruption time with error {e}. Continuing as if the shot did not disrupt.")
               
        mds_conn = self.process_mds_conn.get_shot_connection(shot_id=shot_id)
        
        mds_conn.add_tree_nickname_funcs(
            tree_nickname_funcs = { 
                "_efit_tree" : self.get_efit_tree_nickname_func(
                    shot_id=shot_id, 
                    mds_conn=mds_conn,
                    disruption_time=disruption_time,
                    shot_settings=shot_settings,
                )
            }
        )
        
        try:
            shot_props = self.setup_shot_props(
                shot_id=shot_id,
                mds_conn = mds_conn,
                database=self.process_database,
                disruption_time=disruption_time,
                shot_settings=shot_settings,
                tokamak=Tokamak.CMOD,
                **kwargs
            )
            return shot_props
        except Exception as e:
            self.logger.info(f"[Shot {shot_id}]: Caught failed to setup shot {shot_id}, cleaning up tree manager.")
            mds_conn.close_all_trees()
            raise e
    
    @classmethod
    def get_efit_tree_nickname_func(cls, shot_id : int, mds_conn : MDSConnection, disruption_time : float, shot_settings : ShotSettings) -> None:
        def efit_tree_nickname_func():
            efit_names_to_test = without_duplicates([
                shot_settings.efit_tree_name,
                "analysis",
                *[f"efit0{i}" for i in range(1, 10)],
                *[f"efit{i}" for i in range(10, 19)],
            ])
            
            for efit_name in efit_names_to_test:
                try:
                    mds_conn.open_tree(efit_name)
                    return efit_name
                except Exception as e:
                    cls.logger.info(f"[Shot {shot_id}]: Failed to open efit tree {efit_name} with error {e}.")
                    continue
            
            raise Exception(f"Failed to find efit tree with name {shot_settings.efit_tree_name} in shot {shot_id}.")
        return efit_tree_nickname_func
    
    
    @classmethod
    def _modify_times_disruption_timebase(cls, shot_props : ShotProps, **kwargs):
        raw_ip, ip_time = shot_props.mds_conn.get_record_data(f"ptdata('ip', {shot_props.shot_id})", tree_name='d3d')
        ip_time = ip_time/1.e3
        baseline = np.mean(raw_ip[0:10])
        ip = raw_ip - baseline
        duration, ip_max = cls._get_end_of_shot(ip, ip_time, 100e3)
        if duration < 0.1 or np.abs(ip_max) < 400.e3:
            raise NotImplementedError()
        times = np.arange(0.100, duration+0.025, 0.025)
        if shot_props.disrupted:
            additional_times = np.arange(
                shot_props.disruption_time-cls.DURATION_BEFORE_DISRUPTION, 
                shot_props.disruption_time + cls.DT_BEFORE_DISRUPTION, 
                cls.DT_BEFORE_DISRUPTION
            )
            times = times[np.where(times < (shot_props.disruption_time - cls.DURATION_BEFORE_DISRUPTION))]
            times = np.concatenate((times, additional_times))
        else:
            ip_start = np.argmax(ip_time <= .1)
            ip_end = np.argmax(raw_ip[ip_start:] <= 100000) + ip_start
            return ip_time[ip_start:ip_end]  # [ms] -> [s]
        shot_props.times = times
        return shot_props
        
    @classmethod
    def _modify_times_flattop_timebase(cls, shot_props : ShotProps, **kwargs):
        try:
            ip_prog, t_ip_prog, = shot_props.mds_conn.get_record_data(f"ptdata('iptipp', {shot_props.shot_id})", tree_name='d3d')
            t_ip_prog = t_ip_prog/1.e3  # [ms] -> [s]
            polarity = np.unique(shot_props.mds_conn.get(f"ptdata('iptdirect', {shot_props.shot_id})", tree_name='d3d').data())
            if len(polarity) > 1:
                cls.logger.info(f"[Shot {shot_props.shot_id}]:Polarity of Ip target is not constant. Using value at first timestep.")
                cls.logger.debug(f"[Shot {shot_props.shot_id}]: Polarity array {polarity}")
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, shot_props.times, 'linear')
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, shot_props.times, 'linear')
        except Exception as e:
            cls.logger.info(
                f"[Shot {shot_props.shot_id}]:Failed to get programmed plasma current parameters")
            cls.logger.debug(
                f"[Shot {shot_props.shot_id}]:{traceback.format_exc()}")
        epsoff, t_epsoff = shot_props.mds_conn.get_record_data(f"ptdata('epsoff', {shot_props.shot_id})", tree_name='d3d')
        t_epsoff =t_epsoff/1.e3 + .001  # [ms] -> [s] # Avoid problem with simultaneity of epsoff being triggered exactly on the last time sample
        epsoff = interp1(t_epsoff, epsoff, shot_props.times, 'linear')
        railed_indices = np.where(np.abs(epsoff) > .5)
        power_supply_railed = np.zeros(len(shot_props.times))
        power_supply_railed[railed_indices] = 1
        indices_flattop = np.where((np.abs(dipprog_dt) <= 2.e3) & (
            np.abs(ip_prog) > 100e3) & (power_supply_railed != 1))
        shot_props.times = shot_props.times[indices_flattop]
        return shot_props
    
    @classmethod  
    def _modify_times_rampup_and_flattop_timebase(cls, shot_props : ShotProps, **kwargs):
        raise NotImplementedError()
    
    @classmethod
    def _get_end_of_shot(cls, signal, signal_time, threshold=1.e5):
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