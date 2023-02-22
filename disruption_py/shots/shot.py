from disruption_py.utils import interp1
import subprocess

import MDSplus
from MDSplus import *

import pandas as pd
import numpy as np
import logging

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']


class Shot:
    logger = logging.getLogger('disruption_py')

    def __init__(self, shot_id, data=None):
        self._shot_id = int(shot_id)
        try:
            commit_hash = subprocess.check_output(["git", "describe", "--always"]).strip()
        except Exception as e:
            commit_hash = 'Unknown'
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

    @classmethod
    def get_signal(id, signal_name):
        

    def get_signal(self, signal_name, remote=False, interpolate=True, interpolation_timebase=None):
        """Get a signal from MDSplus.

        Parameters
        ----------
        signal_name : str
            Name of the signal in MDSplus.
        signal_getter : MDSplus.Connection, optional
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
        if signal_getter is None:
            signal_getter = self.conn
            signal = signal_getter.get(signal_name).data()
            orig_timebase = signal_getter.get(
                f"dim_of({signal_name})").data()/1.e3  # [ms] -> [s]
        elif isinstance(MDSplus.Tree, signal_getter):
            signal_record = signal_getter.getNode(signal_name).getData()
            signal = signal_record.data()
            orig_timebase = signal_record.dim_of(0)
        else:
            raise TypeError(
                "signal_getter must be either MDSplus.Connection or MDSplus.Tree")
        if len(orig_timebase) < 2:
            raise ValueError(
                f"Timebase for {signal_name} is too short ({len(orig_timebase)})")
        if interpolate:
            if interpolation_timebase is None:
                interpolation_timebase = self.data['time']
            signal = interp1(orig_timebase, signal, interpolation_timebase)
        return signal, orig_timebase
