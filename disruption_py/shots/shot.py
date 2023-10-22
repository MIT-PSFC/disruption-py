from disruption_py.utils import interp1
from typing import Set
import subprocess
import os
import time
import traceback
import concurrent 
from concurrent.futures import ThreadPoolExecutor

import MDSplus
from MDSplus import *

from disruption_py.mdsplus_integration.tree_manager import TreeManager
from disruption_py.method_caching import MethodOptimizer, CachedMethod

import pandas as pd
import numpy as np
import logging

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']


class Shot:
    """
    Class for a single shot.

    Parameters
    ----------
    shot_id : int
        Shot number.
    data : pandas.DataFrame, optional
        Data for the shot. If not provided, an empty DataFrame will be created.

    Attributes
    ----------
    data : pandas.DataFrame
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

    def __init__(self, shot_id, data=None, **kwargs):
        self._shot_id = int(shot_id)
        self._tree_manager = TreeManager(shot_id)
        self.multiprocessing = kwargs.get('multiprocessing', False)
        print("We are multiprocessing!" if self.multiprocessing else "We are not multiprocessing")
        if self.logger.level == logging.NOTSET:
            self.logger.setLevel(logging.INFO)
        assert self.logger.level != logging.NOTSET, "Logger level is NOTSET"
        if not self.logger.hasHandlers():
            self.logger.addHandler(logging.StreamHandler())
            
        try:
            commit_hash = subprocess.check_output(
                ["git", "describe", "--always"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT).strip()
        except Exception as e:
            # self.logger.warning("Git commit not found")
            commit_hash = 'Unknown'
            
        assert self.logger.hasHandlers(), "Logger has no handlers"
        self._metadata = {
            'labels': {},
            'commit_hash': commit_hash,
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  # TODO: Fix
        }
        self.data = data
        if data is None:
            self.data = pd.DataFrame()

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

    def populate_method(self, cached_method : CachedMethod, method_optimizer : MethodOptimizer, start_time):
        method = getattr(self, cached_method.name)
        result = None
        if callable(method) and hasattr(method, 'populate'):
            self.logger.info(
                f"[Shot {self._shot_id}]:Populating {cached_method.name}")
            # self._tree_manager.cleanup_not_needed()
            try:
                result = method()
            except Exception as e:
                self.logger.warning(
                    f"[Shot {self._shot_id}]:Failed to populate {cached_method.name} with error {e}")
                self.logger.debug(f"{traceback.format_exc()}")
        elif callable(method) and hasattr(method, 'cached'):
            self.logger.info(
                f"[Shot {self._shot_id}]:Caching {cached_method.name}")
            # self._tree_manager.cleanup_not_needed()
            try:
                method()
            except Exception as e:
                self.logger.warning(
                    f"[Shot {self._shot_id}]:Failed to cache {cached_method.name} with error {e}")
                self.logger.debug(f"{traceback.format_exc()}")
        else:
            self.logger.warning(
                f"[Shot {self._shot_id}]:Method {cached_method.name} is not callable or does not have a `populate` attribute set to True")
            return None
        
        print(f"[Shot {self._shot_id}]:Completed {cached_method.name}, time_elapsed: {time.time() - start_time}")
        return result
    def populate_methods(self, method_names):
        """Populate the shot object with data from MDSplus.

        Parameters
        ----------
        method_names : list of str
            List of methods to populate. Each method must be a method of the Shot class
            and must have a `populate` attribute set to True.
        """
        local_data = []
        for method_name in method_names:
            local_data.append(self.populate_method(method_name))
        self.data = pd.concat([self.data, *local_data], axis=1)

    def _init_populate(self, already_populated, methods, tags):
        """
        Internal method to populate the disruption parameters of a shot object. 
        This method is called by the constructor and should not be called directly. It loops through all methods of the Shot class and calls the ones that have a `populate` attribute set to True and satisfy the tags and methods arguments.
        """
        
        # If the shot object was already passed data in the constructor, use that data. Otherwise, create an empty dataframe.
        if not already_populated:
            if self.data is None:
                self.data = pd.DataFrame()
            self.data['time'] = self._times
            self.data['shot'] = self._shot_id
        else:
            self.logger.info(f"[Shot {self._shot_id}]:Already populated")
            return

        # If tags or methods are not lists, make them lists.
        if tags is not None and not isinstance(tags, list):
            populate_tags = [populate_tags]
        if methods is not None and not isinstance(methods, list):
            populate_methods = [populate_methods]
        parameters = []

        # Loop through each attribute and find methods that should populate the shot object.
        methods_to_evaluate : list[CachedMethod] = []
        all_cached_methods : list[CachedMethod] = []
        for method_name in dir(self):
            method = getattr(self, method_name)
            if callable(method) and hasattr(method, 'cached'):
                param_method = CachedMethod(name=method_name, method=method)
                all_cached_methods.append(param_method)
            if callable(method) and hasattr(method, 'populate'):
                param_method = param_method or CachedMethod(name=method_name, method=method)
                # If method does not have tag included and name included then skip
                if tags is not None and bool(set(method.tags).intersection(tags)):
                    methods_to_evaluate.append(param_method)
                    continue
                if methods is not None and method_name in methods:
                    methods_to_evaluate.append(param_method)
                    continue
                self.logger.info(
                        f"[Shot {self._shot_id}]:Skipping {method_name}")
        method_optimizer : MethodOptimizer = MethodOptimizer(self._tree_manager, methods_to_evaluate, all_cached_methods)
        
        if self.multiprocessing:
            # TODO: using method optimizer
            futures = set()
            future_method_names = {}
            def future_for_next(next_method):
                new_future = executor.submit(self.populate_method, next_method, method_optimizer, start_time)
                futures.add(new_future)
                future_method_names[new_future] = next_method.name
            
            start_time = time.time()
            available_methods_runner = method_optimizer.get_async_available_methods_runner(future_for_next)
            with ThreadPoolExecutor(max_workers=1) as executor: 
                available_methods_runner()
                while futures:
                    done, futures = concurrent.futures.wait(futures, return_when='FIRST_COMPLETED')
                    for future in done:
                        try:
                            parameter_df = future.result()
                            parameters.append(parameter_df)
                            method_optimizer.method_complete(future_method_names[future])
                            self._tree_manager.cleanup_not_needed(method_optimizer.can_tree_be_closed)
                        except Exception as e:
                            self.logger.warning(
                                f"[Shot {self._shot_id}]:Failed to populate {method_name} with future error {e}")
                            self.logger.debug(
                                f"[Shot {self._shot_id}: {traceback.format_exc()}")
                    available_methods_runner()
                    
        else:
            parameters = []
            start_time = time.time()
            method_optimizer.run_methods_sync(
                lambda next_method: parameters.append(self.populate_method(next_method, method_optimizer, start_time)),
                self._tree_manager.get_open_tree_names
            )
        parameters = [
            parameter for parameter in parameters if parameter is not None]
        # TODO: This is a hack to get around the fact that some methods return
        #       multiple parameters. This should be fixed in the future.
        local_data = pd.concat(parameters + [self.data], axis=1)
        local_data = local_data.loc[:, ~local_data.columns.duplicated()]
        self.data = local_data
