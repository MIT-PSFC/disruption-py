from disruption_py.utils import interp1
import subprocess

import MDSplus
from MDSplus import *

import pandas as pd

DEFAULT_SHOT_COLUMNS = ['time', 'shot', 'time_until_disrupt', 'ip']


class Shot:
    def __init__(self, shot_id, data_columns, data=None):
        self._shot_id = shot_id
        self._metadata = {
            'labels': {},
            'commit_hash': subprocess.check_output(["git", "describe", "--always"]).strip(),
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  # TODO: Fix
        }
        self.data = data
        if data is None:
            self.data = pd.DataFrame()

    def get_signal(self, signal_name, signal_getter, interpolate=True, interpolation_timebase=None):
        if isinstance(MDSplus.Connection, signal_getter):
            signal = signal_getter.get(signal_name).data()
            orig_timebase = signal_getter.get(
                f"dim_of({signal_name})").data()/1e3  # [ms] -> [s]
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
