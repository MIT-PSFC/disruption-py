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
            commit_hash = subprocess.check_output(
                ["git", "describe", "--always"]).strip()
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
    def get_signal(id, signal_name, conn, interpolate=True, interpolation_timebase=None):
        if isinstance(MDSplus.Tree, conn):
            signal_record = conn.getNode(signal_name).getData()
            signal = signal_record.data()
            orig_timebase = signal_record.dim_of(0)
        elif isinstance(MDSplus.Connection, conn):
            signal = conn.get(signal_name).data()
            orig_timebase = conn.get(
                f"dim_of({signal_name})").data()/1.e3  # [ms] -> [s]
        else:
            raise TypeError(
                "conn must be either MDSplus.Connection or MDSplus.Tree")
        if len(orig_timebase) < 2:
            raise ValueError(
                f"Timebase for {signal_name} is too short ({len(orig_timebase)})")
        if interpolate:
            if interpolation_timebase is None:
                raise ValueError(
                    "interpolation_timebase must be provided if interpolate is True")
            signal = interp1(orig_timebase, signal, interpolation_timebase)
        return signal, orig_timebase

    def _get_signal(self, signal_name, signal_conn=None, interpolate=True, interpolation_timebase=None):
        """Get a signal from MDSplus.

        Parameters
        ----------
        signal_name : str
            Name of the signal in MDSplus.
        signal_conn : MDSplus.Connection, optional
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
        if signal_conn is None:
            signal_conn = self.conn
        if interpolation_timebase is None and interpolate:
            interpolation_timebase = self.data['time']
        return type(self).get_signal(self._shot_id, signal_name, signal_conn, interpolate, interpolation_timebase)
