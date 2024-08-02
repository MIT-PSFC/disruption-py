#!/usr/bin/env python3

import traceback
from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Dict, Union

import numpy as np
import pandas as pd

from disruption_py.config import config
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.io.mds import MDSConnection
from disruption_py.io.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class TimeSettingParams:
    """Params passed by disruption_py to _get_times() method.

    Attributes
    ----------
    shot_id : int
        The shot id for the timebase being created
    mds_conn : MDSConnection
        The connection to MDSplus, that can be used to retrieve MDSplus data.
    cache_data : pd.DataFrame
        Pre-filled data given to disruption_py.
    database : ShotDatabase
        Database object with connection to sql database.
    disruption_time : float
        The time when the shot disrupted or None if no disruption occured.
    tokamak : Tokamak
        The tokamak using the time setting.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    shot_id: int
    mds_conn: MDSConnection
    cache_data: pd.DataFrame
    database: ShotDatabase
    disruption_time: float
    tokamak: Tokamak
    logger: Logger

    @property
    def disrupted(self) -> bool:
        """Returns true if the shot disrupted."""
        return self.disruption_time is not None


TimeSettingType = Union["TimeSetting", str, Dict[Tokamak, "TimeSettingType"]]


class TimeSetting(ABC):
    """TimeSetting abstract class that should be inherited by all time setting classes."""

    def get_times(self, params: TimeSettingParams) -> np.ndarray:
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_times(params)

    @abstractmethod
    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """Abstract method implemented by subclasses to get timebase as list.
        The timebase can be set to be automatically restricted to a subdomain of the
        provided times via the domain_setting argument in the RetrievalSettings object.

        Parameters
        ----------
        params : TimeSettingParams
            Params that can be used to determine and retrieve the timebase.

        Returns
        -------
        np.ndarray
            Numpy array containing times in the timebase.
        """
        pass


class TimeSettingDict(TimeSetting):
    """
    Utility class that is automatically used when a dicationary is passed as the `time_setting` parameter in `RetrievalSettings`.

    Parameters
    ----------
    time_setting_dict : dict[Tokamak, TimeSettingType]
        A dictionary mapping tokamak type strings to the desired time setting for that tokamak.  E.g. `{'cmod': 'efit'}`.
        Any other option passable to the `time_setting` parameter in `RetrievalSettings` may be used.
    """

    def __init__(self, time_setting_dict: Dict[Tokamak, TimeSettingType]):
        resolved_time_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_time_setting(
                individual_setting
            )
            for tokamak, individual_setting in time_setting_dict.items()
        }
        self.resolved_time_setting_dict = resolved_time_setting_dict

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        chosen_setting = self.resolved_time_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_times(params)
        else:
            params.logger.warning(f"No time setting for tokamak {params.tokamak}")
            return None


class ListTimeSetting(TimeSetting):
    """
    A list of times to use as the timebase.

    A utility class that is automatically used when a list, numpy array, or pandas series is
    passed as the `time_setting` parameter in `RetrievalSettings`.
    """

    def __init__(self, times):
        self.times = times

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        return self.times


class CacheTimeSetting(TimeSetting):
    """
    Time setting for using the cached data timebase.

    If no cached data exists for the shot, then the fallback_time_setting is used.
    """

    def __init__(self, fallback_time_setting: TimeSettingType) -> None:
        self.fallback_time_setting = resolve_time_setting(fallback_time_setting)

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        if params.cache_data is not None:
            # set timebase to be the timebase of cached data
            try:
                times = params.cache_data["time"].to_numpy()
                # Check if the timebase is in ms instead of s
                if times[-1] > config(params.tokamak).MAX_SHOT_TIME:
                    times /= 1000  # [ms] -> [s]
                return times
            except KeyError as e:
                params.logger.warning(
                    f"[Shot {params.shot_id}]: Shot constructor was passed data but no timebase."
                )
                params.logger.debug(f"[Shot {params.shot_id}]:{traceback.format_exc()}")
        else:
            return self.fallback_time_setting.get_times(params)


class EfitTimeSetting(TimeSetting):
    """Time setting for using the EFIT timebase."""

    def __init__(self):
        self.tokamak_overrides = {Tokamak.CMOD: self.cmod_times}

    def cmod_times(self, params: TimeSettingParams):
        efit_tree_name = params.mds_conn.get_tree_name_of_nickname("_efit_tree")
        if efit_tree_name == "analysis":
            try:
                return params.mds_conn.get_data(
                    r"\analysis::efit_aeqdsk:time",
                    tree_name="_efit_tree",
                    astype="float64",
                )
            except Exception as e:
                return params.mds_conn.get_data(
                    r"\analysis::efit:results:a_eqdsk:time",
                    tree_name="_efit_tree",
                    astype="float64",
                )
        else:
            return params.mds_conn.get_data(
                rf"\{efit_tree_name}::top.results.a_eqdsk:time",
                tree_name="_efit_tree",
                astype="float64",
            )

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        raise ValueError("EFIT timebase setting not implemented")


class DisruptionTimeSetting(TimeSetting):
    """
    Time setting for using the disruption timebase.

    The disruption timebase is the timebase of the shot that was disrupted.
    """

    # Disruption Variables
    DT_BEFORE_DISRUPTION = 0.002
    DURATION_BEFORE_DISRUPTION = 0.10

    def __init__(self, minimum_ip=400.0e3, minimum_duration=0.100):
        self.tokamak_overrides = {Tokamak.D3D: self.d3d_times}
        self.minimum_ip = minimum_ip
        self.minimum_duration = minimum_duration

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        raise ValueError("Disruption time setting not implemented")

    def d3d_times(self, params: TimeSettingParams):
        raw_ip, ip_time = params.mds_conn.get_data_with_dims(
            f"ptdata('ip', {params.shot_id})", tree_name="d3d"
        )
        ip_time = ip_time / 1.0e3
        baseline = np.mean(raw_ip[0:10])
        ip = raw_ip - baseline
        duration, ip_max = self._get_end_of_shot(ip, ip_time, 100e3)
        if duration < self.minimum_duration or np.abs(ip_max) < self.minimum_ip:
            raise NotImplementedError()

        times = np.arange(0.100, duration + config(params.tokamak).TIME_CONST, 0.025)
        if params.disrupted:
            additional_times = np.arange(
                params.disruption_time - self.DURATION_BEFORE_DISRUPTION,
                params.disruption_time + config(params.tokamak).TIME_CONST,
                self.DT_BEFORE_DISRUPTION,
            )
            times = times[
                np.where(
                    times
                    < (
                        params.disruption_time
                        - self.DURATION_BEFORE_DISRUPTION
                        - config(params.tokamak).TIME_CONST
                    )
                )
            ]
            times = np.concatenate((times, additional_times))
        # else:
        #     ip_start = np.argmax(ip_time <= .1)
        #     ip_end = np.argmax(raw_ip[ip_start:] <= 100000) + ip_start
        #     times = ip_time[ip_start:ip_end]  # [ms] -> [s]
        return times

    @classmethod
    def _get_end_of_shot(cls, signal, signal_time, threshold=1.0e5):
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
        finite_indices = np.where((signal_time >= 0.0) & (np.abs(signal) > threshold))
        if len(finite_indices) == 0:
            return duration, signal_max
        else:
            dt = np.diff(signal_time)
            duration = np.sum(dt[finite_indices[:-1]])
            if duration < 0.1:  # Assuming < 100 ms is not a bona fide plasma
                duration = 0
                return duration, signal_max
        polarity = np.sign(
            np.trapz(signal[finite_indices], signal_time[finite_indices])
        )
        polarized_signal = polarity * signal
        valid_indices = np.where((polarized_signal >= threshold) & (signal_time > 0.0))
        duration = signal_time[np.max(valid_indices)]
        if len(valid_indices) == signal_time.size:
            duration = -duration
        signal_max = np.max(polarized_signal) * polarity
        return duration, signal_max


class IpTimeSetting(TimeSetting):
    def __init__(self):
        self.tokamak_overrides = {Tokamak.D3D: self.d3d_times}

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        raise ValueError("Ip time setting not implemented")

    def d3d_times(self, params: TimeSettingParams):
        (ip_time,) = params.mds_conn.get_dims(
            f"ptdata('ip', {params.shot_id})", tree_name="d3d"
        )
        return ip_time
class PCSTimeSetting(TimeSetting):
    """Plasma Control System Time Setting, runs at 1 kHz"""
    def __init__(self):
        self.tokamak_overrides = {Tokamak.CMOD: self.cmod_times}

    def cmod_times(self, params: TimeSettingParams):
        return np.array(np.arange(-0.1, 2, 0.001)) # TODO(ZanderKeith) Just hardcoded start and stop time for now

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        return self.tokamak_overrides[params.tokamak](params)
    
class OfflineTimeSetting(TimeSetting):
    """Plasma Control System Time Setting, runs at 1 kHz"""
    def __init__(self):
        self.tokamak_overrides = {Tokamak.CMOD: self.cmod_times}

    def cmod_times(self, params: TimeSettingParams):
        return np.array(np.arange(0.6, 1.6, 0.0001)) # TODO(ZanderKeith) Just hardcoded start and stop time for now
        
    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        return self.tokamak_overrides[params.tokamak](params)

class SignalTimeSetting(TimeSetting):
    def __init__(self, tree_name: str, signal_path: str):
        self.tree_name = tree_name
        self.signal_path = signal_path

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        try:
            (signal_time,) = params.mds_conn.get_dims(
                self.signal_path, tree_name=self.tree_name, astype="float64"
            )
            return signal_time
        except Exception as e:
            params.logger.error(
                f"Failed to set up timebase for signal {self.signal_path}"
            )
            raise Exception(e)


# --8<-- [start:time_setting_dict]
_time_setting_mappings: Dict[str, TimeSetting] = {
    "efit": EfitTimeSetting(),
    "disruption": DisruptionTimeSetting(),
    "disruption_warning": {
        Tokamak.CMOD: EfitTimeSetting(),
        Tokamak.D3D: DisruptionTimeSetting(),
    },
    "ip": IpTimeSetting(),
    "pcs": PCSTimeSetting(),
    "offline": OfflineTimeSetting(),
}
# --8<-- [end:time_setting_dict]


def resolve_time_setting(
    time_setting: TimeSettingType,
) -> TimeSetting:
    if isinstance(time_setting, TimeSetting):
        return time_setting

    if isinstance(time_setting, str):
        time_setting_object = _time_setting_mappings.get(time_setting, None)
        if time_setting_object is not None:
            return time_setting_object

    if isinstance(time_setting, np.ndarray) or isinstance(time_setting, list):
        return ListTimeSetting(time_setting)

    if isinstance(time_setting, pd.Series):
        return ListTimeSetting(time_setting.to_numpy())

    if isinstance(time_setting, dict):
        return TimeSettingDict(time_setting)

    raise ValueError("Invalid time setting")
