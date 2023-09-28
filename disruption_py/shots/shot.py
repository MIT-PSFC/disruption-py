from disruption_py.utils import interp1
import subprocess
import traceback
import concurrent 
from concurrent.futures import ThreadPoolExecutor

import MDSplus
from MDSplus import *

from disruption_py.mdsplus_integration.tree_manager import TreeManager

import pandas as pd
import numpy as np
import logging

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']


def parameter_method(tags=["all"]):
    """
    Tags a function as a parameter method and instantiates its cache. Parameter methods are functions that 
    calculate disruption parameters from the data in the shot.  They are called by the Shot object when
    it is instantiated. The cache is used to store the results of the parameter method so that it is only
    calculated once per shot for a given timebase.
    """
    # TODO: Figure out how to hash _times so that we can use the cache for different timebases
    def tag_wrapper(func):
        def wrapper(self, *args, **kwargs):
            # Create the cache if it doesn't exist
            if not hasattr(self, '_cached_result'):
                self._cached_result = {}
            cache_key = func.__name__ + str(len(self._times))
            if cache_key in self._cached_result:
                return self._cached_result[cache_key]
            else:
                result = func(self, *args, **kwargs)
                self._cached_result[cache_key] = result
                return result

        wrapper.populate = True
        wrapper.tags = tags
        return wrapper
    return tag_wrapper


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
        try:
            commit_hash = subprocess.check_output(
                ["git", "describe", "--always"]).strip()
        except Exception as e:
            commit_hash = 'Unknown'
        if self.logger.level == logging.NOTSET:
            self.logger.setLevel(logging.INFO)
        assert self.logger.level != logging.NOTSET, "Logger level is NOTSET"
        if not self.logger.hasHandlers():
            self.logger.addHandler(logging.StreamHandler())
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

    def populate_method(self, method_name):
        method = getattr(self, method_name)
        if callable(method) and hasattr(method, 'populate'):
            self.logger.info(
                f"[Shot {self._shot_id}]:Populating {method_name}")
            try:
                return method()
            except Exception as e:
                self.logger.warning(
                    f"[Shot {self._shot_id}]:Failed to populate {method_name}")
                self.logger.debug(f"{traceback.format_exc()}")
        else:
            self.logger.warning(
                f"[Shot {self._shot_id}]:Method {method_name} is not callable or does not have a `populate` attribute set to True")
        return None
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

    def populate_tag(self, tag):
        local_data = []
        for method_name in dir(self):
            method = getattr(self, method_name)
            if callable(method) and hasattr(method, 'populate'):
                if bool(set(method.tag).intersection(tag)):
                    print(f"[Shot {self._shot_id}]:Skipping {method_name}")
                    continue
                try:
                    local_data.append(method())
                except Exception as e:
                    self.logger.warning(
                        f"[Shot {self._shot_id}]:Failed to populate {method_name}")
                    self.logger.debug(f"{traceback.format_exc()}")
        self.data = pd.concat([self.data, local_data], axis=1)

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
        method_names = []
        for method_name in dir(self):
            method = getattr(self, method_name)
            if callable(method) and hasattr(method, 'populate'):
                # If method does not have tag included and name included then skip
                if tags is not None and bool(set(method.tags).intersection(tags)):
                    method_names.append(method_name)
                    continue
                if methods is not None and method_name in methods:
                    method_names.append(method_name)
                    continue
                self.logger.info(
                        f"[Shot {self._shot_id}]:Skipping {method_name}")
        if self.multiprocessing:
            with ThreadPoolExecutor(max_workers=8) as executor:
                futures = [executor.submit(
                    self.populate_method, method_name) for method_name in method_names]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        parameter_df = future.result()
                        parameters.append(parameter_df)
                    except Exception as e:
                        self.logger.warning(
                            f"[Shot {self._shot_id}]:Failed to populate {method_name}")
                        self.logger.debug(
                            f"[Shot {self._shot_id}: {traceback.format_exc()}")
        else:
            parameters = [self.populate_method(
                method_name) for method_name in method_names]
        parameters = [
            parameter for parameter in parameters if parameter is not None]
        # TODO: This is a hack to get around the fact that some methods return
        #       multiple parameters. This should be fixed in the future.
        local_data = pd.concat(parameters + [self.data], axis=1)
        local_data = local_data.loc[:, ~local_data.columns.duplicated()]
        self.data = local_data
