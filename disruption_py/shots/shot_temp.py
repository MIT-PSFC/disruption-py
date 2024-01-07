from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Tuple
from disruption_py.settings.enum_options import InterpolationMethod
from disruption_py.utils.command_utils import get_commit_hash
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.mappings.tokamak import Tokamak
import subprocess
import traceback


import MDSplus
from MDSplus import *

from disruption_py.mdsplus_integration.tree_manager import TreeManager, EnvModifications
from disruption_py.settings.enum_options import SignalDomain
from disruption_py.settings.set_times_request import SetTimesRequest, SetTimesRequestParams
from disruption_py.utils.constants import TIME_CONST

import pandas as pd
import numpy as np
import logging

from disruption_py.utils.utils import without_duplicates

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']
MAX_SHOT_TIME = 7.0  # [s]

@dataclass
class ShotSetupParams:
    shot_id : int
    tokamak: Tokamak
    num_threads_per_shot : int
    override_exising_data : bool
    set_times_request : SetTimesRequest
    signal_domain : SignalDomain
    existing_data : pd.DataFrame
    disruption_time : float
    tree_nicknames : Dict[str, Tuple[List[str], List[EnvModifications]]]

class Shot(ABC):
    logger = logging.getLogger('disruption_py')

    def __init__(
        self, 
        setup_params : ShotSetupParams=None,
        **kwargs
    ):
        self._shot_id = setup_params.shot_id
        self._tokamak = setup_params.tokamak
        self._num_threads_per_shot = setup_params.num_threads_per_shot
        self._disruption_time = setup_params.disruption_time
        self._tree_manager = TreeManager(setup_params.shot_id)
        self._initial_existing_data = setup_params.existing_data
            
        if self._num_threads_per_shot > 1:
            self.logger.info("Intra-shot multithreading enabled")
            
        self._metadata = {
            'labels': {},
            'commit_hash': get_commit_hash(),
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  # TODO: Fix
        }
        
        self._setup_nicknames(setup_params.tree_nicknames)

        # must call this method after nicknamed trees setup
        self._init_timebase(setup_params, setup_params.existing_data)
        
        # add existing data
        self._init_with_data(setup_params.existing_data)
    
    @abstractmethod
    def set_flattop_timebase(self):
        pass
    
    @abstractmethod
    def set_rampup_and_flattop_timebase(self):
        pass
    
    def cleanup(self):
        """
        Remove references to mdsplus resources to prevent data leaks.
        """
        self._tree_manager.cleanup()
        self._times = None
        if hasattr(self, '_cached_result'):
            self._cached_result.clear()

    @property
    def times(self):
        return self._times
    
    @property
    def shot_id(self):
        return self._shot_id
    
    @property
    def commit_hash(self):
        return self._metadata.get('commit_hash', 'Unknown')
    
    @property
    def tree_manager(self):
        return self._tree_manager
    
    @property
    def disrupted(self):
        return self._disruption_time is not None
    
    @property
    def disruption_time(self):
        return self._disruption_time
    
    @property
    def initial_existing_data(self):
        return self._initial_existing_data
    
    @property
    def metadata(self):
        return self._metadata
    
    def __getitem__(self, key):
        return self._metadata if key == 'metadata' else self.data[key]
    
    def _init_timebase(self, setup_params: ShotSetupParams, existing_data):
        """
        Initialize the timebase of the shot.
        """
        if existing_data is not None and setup_params.override_exising_data is False:
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
            request_params = SetTimesRequestParams(tree_manager=self._tree_manager, tokamak=self._tokamak, logger=self.logger, disruption_time=self._disruption_time)
            self._times = setup_params.set_times_request.get_times(request_params)
        self.interpolation_method : InterpolationMethod  = interp1 # TODO: fix
        
        if setup_params.signal_domain is SignalDomain.FLATTOP:
            self.set_flattop_timebase()
        elif setup_params.signal_domain is SignalDomain.RAMP_UP_AND_FLATTOP:
            self.set_rampup_and_flattop_timebase()
            
    def _setup_nicknames(self, tree_nicknames : Dict[str, Tuple[List[str], List[EnvModifications]]]):
        # nickname efit tree
        for nickname, (tree_names, env_modifications) in tree_nicknames.items():
            efit_names_to_test = without_duplicates(tree_names)
            self._tree_manager.nickname(nickname, efit_names_to_test, env_modifications)
    
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

    