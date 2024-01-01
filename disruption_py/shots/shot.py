from abc import ABC, abstractmethod
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.mappings.tokamak import Tokamak
import subprocess
import traceback


import MDSplus
from MDSplus import *

from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.settings.shot_settings import ShotSettings, InterpolationMethod, SignalDomain
from disruption_py.settings.set_times_request import SetTimesRequestParams
from disruption_py.utils.constants import TIME_CONST

import pandas as pd
import numpy as np
import logging

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']
MAX_SHOT_TIME = 7.0  # [s]


class Shot(ABC):
    """
    Abstract class for a single shot.

    Parameters
    ----------
    shot_id : int
        Shot number.
    existing_data : pandas.DataFrame, optional
        Data for the shot. If not provided, an empty DataFrame will be created.

    Attributes
    ----------
    existing_data : pandas.DataFrame
        Data for the shot.
    conn : MDSplus.Connection
        MDSplus connection to the shot.
    tree : MDSplus.Tree
        MDSplus tree for the shot.
    logger : logging.Logger
        Logger for the shot.
    _shot_id : int
        Shot number.
    _metadata : dict
        Metadata for the shot.
    """
    # TODO: Add [Shot {self._shot_id}]: to logger format by default
    logger = logging.getLogger('disruption_py')

    def __init__(
        self, 
        shot_id, 
        tokamak: Tokamak,
        disruption_time=None,
        shot_settings : ShotSettings=None,
        **kwargs
    ):
        self._shot_id = int(shot_id)
        self.tokamak = tokamak
        self.num_threads_per_shot = shot_settings.num_threads_per_shot
        self.disruption_time = disruption_time
        self.disrupted = self.disruption_time is not None
        self._tree_manager = TreeManager(shot_id)        
            
        if self.num_threads_per_shot > 1:
            self.logger.info("Multithreading enabled")
        
        # setup commit hash
        try:
            commit_hash = subprocess.check_output(
                ["git", "describe", "--always"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT).strip()
        except Exception as e:
            # self.logger.warning("Git commit not found")
            commit_hash = 'Unknown'
            
        self._metadata = {
            'labels': {},
            'commit_hash': commit_hash,
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  # TODO: Fix
        }
        
    @abstractmethod
    def setup_nicknames(self):
        pass
    
    @abstractmethod
    def set_flattop_timebase(self):
        pass
    
    @abstractmethod
    def set_rampup_and_flattop_timebase(self):
        pass
    
    @abstractmethod
    def cleanup(self):
        pass

    def get_times(self):
        return self._times
    
    def get_shot_id(self):
        return self._shot_id
    
    def get_commit_hash(self):
        return self._metadata.get('commit_hash', 'Unknown')
    
    def get_tree_manager(self):
        return self._tree_manager
    
    def get_disruption_time(self):
        return self.disruption_time
    
    def _init_timebase(self, shot_settings: ShotSettings, existing_data):
        """
        Initialize the timebase of the shot.
        """
        if existing_data is not None and shot_settings.override_exising_data is False:
            # set timebase to be the timebase of existing data
            try:
                self._times = existing_data['time'].to_numpy()
                # Check if the timebase is in ms instead of s
                if self._times[-1] > MAX_SHOT_TIME:
                    self._times /= 1000  # [ms] -> [s]
            except KeyError as e:
                self.logger.warning(
                    f"[Shot {self._shot_id}]: Shot constructor was passed data but no timebase.")
                self.logger.debug(
                    f"[Shot {self._shot_id}]:{traceback.format_exc()}")
        else:
            request_params = SetTimesRequestParams(tree_manager=self._tree_manager, tokamak=self.tokamak, logger=self.logger, disruption_time=self.disruption_time)
            self._times = shot_settings.set_times_request.get_times(request_params)
        self.interpolation_method : InterpolationMethod  = interp1 # TODO: fix
        
        if shot_settings.signal_domain is SignalDomain.FLATTOP:
            self.set_flattop_timebase()
        elif shot_settings.signal_domain is SignalDomain.RAMP_UP_AND_FLATTOP:
            self.set_rampup_and_flattop_timebase()
    
    def _init_with_data (self, existing_data : pd.DataFrame):
        '''
        Intialize the shot with data, if existing data matches the shot timebase.
        '''
        if existing_data is not None:
            time_df = pd.DataFrame(self._times, columns=['time'])
            flagged_existing_data = existing_data.assign(merge_success_flag=1)
            timed_existing_data = pd.merge_asof(time_df, flagged_existing_data, on='time', direction='nearest', tolerance=TIME_CONST)
            if timed_existing_data['merge_success_flag'].isna().any():
                existing_data = None
            else:
                existing_data = timed_existing_data.drop(columns=['merge_success_flag'])
        
        self.initialized_with_data = existing_data is not None
        self.data = existing_data

    @staticmethod
    def get_signal(signal, conn, interpolate=True, interpolation_timebase=None):
        """Get a signal from MDSplus.

        Parameters
        ----------
        signal : str
            Name of the signal in MDSplus.
        conn : MDSplus.Connection, optional
            MDSplus connection to get the signal from. If not provided, the default
            connection will be used.
        interpolate : bool, optional
            Whether to interpolate the signal onto the timebase of the experiment.
            If True, the signal will be interpolated from the timebase of the signal
            in MDSplus to the timebase of the experiment.
        interpolation_timebase : array_like, optional
            Timebase to interpolate the signal to. If not provided, the current timebase of
            the Shot object will be used.

        Returns
        -------
        signal : array_like
            Signal from MDSplus.
        orig_timebase : array_like
            Timebase of the signal in MDSplus.

        """
        if isinstance(conn, MDSplus.Tree):
            signal_record = conn.getNode(signal).getData()
            signal_data = signal_record.data()
            signal_data = signal_record.data()
            orig_timebase = signal_record.dim_of(0)
        elif isinstance(conn, MDSplus.Connection):
            signal_data = conn.get(signal).data()
            signal_data = conn.get(signal).data()
            orig_timebase = conn.get(
                f"dim_of({signal})").data()/1.e3  # [ms] -> [s]
        else:
            raise TypeError(
                "conn must be either MDSplus.Connection or MDSplus.Tree")
        if len(orig_timebase) < 2:
            raise ValueError(
                f"Timebase for {signal} is too short ({len(orig_timebase)})")
        if interpolate:
            if interpolation_timebase is None:
                raise ValueError(
                    "interpolation_timebase must be provided if interpolate is True")
            signal_data = interp1(
                orig_timebase, signal_data, interpolation_timebase)
        return signal_data, orig_timebase

    def _get_signal(self, signal, conn=None, interpolate=False, interpolation_timebase=None):
        if conn is None:
            conn = self.conn
        if interpolation_timebase is None and interpolate:
            interpolation_timebase = self._times
        return type(self).get_signal(signal, conn, interpolate, interpolation_timebase)

    @staticmethod
    def get_signals(signals, conn, interpolate=True, interpolation_timebase=None):
        """ Get multiple signals from MDSplus. Signal names are executed in order meaning that while this method expects all signals to be on the same tree you can grab signals from multiple trees by adding open tree TDI expressions as signals 

        Parameters
        ----------
        signals : list of str
            Names of the signals in MDSplus.
        conn : MDSplus.Connection
            MDSplus connection to get the signals from.
        interpolate : bool, optional
            Whether to interpolate the signals onto the timebase of the experiment.
            If True, the signals will be interpolated from the timebase of the signals
            in MDSplus to the timebase of the experiment.
        interpolation_timebase : array_like, optional
            Timebase to interpolate the signals to. If not provided, the current timebase of
            the Shot object will be used.

        Returns
        -------
        signals : list of array_like
            Signals from MDSplus.
        orig_timebases : list of array_like
            Timebases of the signals in MDSplus.
        """
        gm = MDSplus.Data.GetMany(conn)
        for path in signals:
            gm.append(path)

        results = gm.get()
        for i, result in enumerate(results):
            if isinstance(result, MDSplus.Data.Error):
                raise ValueError(
                    f"Error getting signal {signals[i]}: {result}")
            if interpolate:
                if interpolation_timebase is None:
                    raise ValueError(
                        "interpolation_timebase must be provided if interpolate is True")
                results[i] = interp1(result.dim_of(
                    0).data(), result.data(), interpolation_timebase)
        return [result.data() for result in results], [result.dim_of(0).data() for result in results]

    def _get_signals(self, signals, conn=None, interpolate=False, interpolation_timebase=None):
        if conn is None:
            conn = self.conn
        if interpolation_timebase is None and interpolate:
            interpolation_timebase = self._times
        return type(self).get_signals(signals, conn, interpolate, interpolation_timebase)
    
    def apply_shot_filter(self, shot_filter):
        self.data = self.data.filter(shot_filter)

    def apply_shot_transform(self, shot_transform):
        self.data = self.data.apply(shot_transform)

    