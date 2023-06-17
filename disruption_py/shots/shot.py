from disruption_py.utils import interp1
import subprocess
import traceback
import traceback

import MDSplus
from MDSplus import *

import pandas as pd
import numpy as np
import logging

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']



def parameter_method(func, tags=["all"]):
    func.populate = True
    func.populate = True
    func.tags = tags
    return func
    return func


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

    logger = logging.getLogger('disruption_py')

    def __init__(self, shot_id, data=None):
        self._shot_id = int(shot_id)
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
            orig_timebase = signal_record.dim_of(0)
        elif isinstance(conn, MDSplus.Connection):
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
            signal_data = interp1(
                orig_timebase, signal_data, interpolation_timebase)
        return signal_data, orig_timebase

    def _get_signal(self, signal, conn=None, interpolate=True, interpolation_timebase=None):
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

    def apply_shot_filter(self, shot_filter):
        self.data = self.data.filter(shot_filter)


    def apply_shot_transform(self, shot_transform):
        self.data = self.data.apply(shot_transform)

    def populate_methods(self, methods_to_populate):
        """Populate the shot object with data from MDSplus.

        Parameters
        ----------
        methods_to_populate : list of str
            List of methods to populate. Each method must be a method of the Shot class
            and must have a `populate` attribute set to True.
        """
        local_data = []
        for method_name in methods_to_populate:
            method = getattr(self, method_name)
            if callable(method) and hasattr(method, 'populate'):
                self.logger.info(f"Populating {method_name}")
                try:
                    local_data.append(method())
                except Exception as e:
                    self.logger.warning(f"Failed to populate {method_name}")
                    self.logger.debug(f"{traceback.format_exc()}")
            else:
                self.logger.warning(
                    f"Method {method_name} is not callable or does not have a `populate` attribute set to True")
        self.data = pd.concat([self.data, *local_data], axis=1)

    def populate_tags(self, tags):
        local_data = []
        for method_name in dir(self):
            method = getattr(self, method_name)
            if callable(method) and hasattr(method, 'populate'):
                if bool(set(method.tags).intersection(tags)):
                    print(f"[Shot {self._shot_id}]:Skipping {method_name}")
                    continue
                try:
                try:
                    local_data.append(method())
                except Exception as e:
                    self.logger.warning(
                        f"[Shot {self._shot_id}]:Failed to populate {method_name}")
                    self.logger.warning(
                        f"[Shot {self._shot_id}]:Failed to populate {method_name}")
                    self.logger.debug(f"{traceback.format_exc()}")
        self.data = pd.concat([self.data, local_data], axis=1)

    def _init_populate(self, already_populated, methods, tags):
    def _init_populate(self, already_populated, methods, tags):
        if not already_populated:
            if self.data is None:
                self.data = pd.DataFrame()
            self.data['time'] = self._times
            self.data['shot'] = self._shot_id
        if tags is not None and not isinstance(tags, list):
            populate_tags = [populate_tags]
        if methods is not None and not isinstance(methods, list):
            populate_methods = [populate_methods]
        parameters = []
        for method_name in dir(self):
            method = getattr(self, method_name)
            if callable(method) and hasattr(method, 'populate'):
                if not bool(set(method.tags).intersection(tags)) or (methods is not None and method_name not in methods):
                    continue
                try:
                try:
                    parameters.append(method())
                except Exception as e:
                    print(e)
                    self.logger.warning(
                        f"[Shot {self._shot_id}]:Failed to populate {method_name}")
                    self.logger.debug(f"{traceback.format_exc()}")
        # TODO: This is a hack to get around the fact that some methods return
        # TODO: This is a hack to get around the fact that some methods return
        #       multiple parameters. This should be fixed in the future.
        local_data = pd.concat(parameters + [self.data], axis=1)
        local_data = local_data.loc[:, ~local_data.columns.duplicated()]
        self.data = local_data
        self.data = local_data
