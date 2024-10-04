#!/usr/bin/env python3

"""
This module defines classes for time settings, used to manage the timebase for 
retrieving data in disruption_py for various tokamaks and shot configurations.
"""

import traceback
from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Dict, Union

import numpy as np
import pandas as pd
from MDSplus import mdsExceptions

from disruption_py.config import config
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.inout.mds import MDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class TimeSettingParams:
    """
    Parameters passed by disruption_py to the _get_times() method.

    Attributes
    ----------
    shot_id : int
        Shot ID for the timebase being created.
    mds_conn : MDSConnection
        Connection to MDSPlus for retrieving MDSPlus data.
    cache_data : pd.DataFrame
        Pre-filled data provided to disruption_py.
    database : ShotDatabase
        Database object with connection to the SQL database.
    disruption_time : float
        Time when the shot disrupted in seconds (or None if no disruption occurred).
    tokamak : Tokamak
        Tokamak for which the time setting is applied.
    logger : Logger
        Logger object used by disruption_py for logging messages.
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
        """
        Check if the shot disrupted.

        Returns
        -------
        bool
            True if the shot was disrupted, False otherwise.
        """
        return self.disruption_time is not None


TimeSettingType = Union["TimeSetting", str, Dict[Tokamak, "TimeSettingType"]]


class TimeSetting(ABC):
    """
    Abstract base class for managing time settings to retrieve the timebase for shots.

    Methods
    -------
    get_times(params)
        Retrieve the timebase as a numpy array.
    """

    def get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Retrieve the timebase using the provided parameters.

        Parameters
        ----------
        params : TimeSettingParams
            The parameters used to determine and retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_times(params)

    @abstractmethod
    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Abstract method for subclasses to implement timebase retrieval.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters used to determine the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """


class TimeSettingDict(TimeSetting):
    """
    Utility class used when a dictionary is passed as the `time_setting` parameter.

    Parameters
    ----------
    time_setting_dict : dict[Tokamak, TimeSettingType]
        Dictionary mapping tokamaks to specific time settings, e.g., `{'cmod': 'efit'}`.
        Any other option passable to the `time_setting` parameter in
        `RetrievalSettings` may be used.
    """

    def __init__(self, time_setting_dict: Dict[Tokamak, TimeSettingType]):
        """
        Initialize with a dictionary mapping tokamaks to time settings.

        Parameters
        ----------
        time_setting_dict : dict
            Dictionary of tokamak to time setting mappings.
        """
        resolved_time_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_time_setting(
                individual_setting
            )
            for tokamak, individual_setting in time_setting_dict.items()
        }
        self.resolved_time_setting_dict = resolved_time_setting_dict

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Retrieve the timebase for the given tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        chosen_setting = self.resolved_time_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_times(params)
        params.logger.warning("No time setting for tokamak %s", params.tokamak)
        return None


class ListTimeSetting(TimeSetting):
    """
    Time setting for using a pre-defined list of times.

    Used when a list, numpy array, or pandas series is passed as the `time_setting`
    parameter in `RetrievalSettings`.
    """

    def __init__(self, times):
        """
        Initialize with a list of times.

        Parameters
        ----------
        times : list, np.ndarray, pd.Series
            List or array of times to use as the timebase.
        """
        self.times = times

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Return the pre-defined list of times.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        return self.times


class CacheTimeSetting(TimeSetting):
    """
    Time setting for using the timebase from cached data.

    If no cached data is available, the fallback time setting is used.
    """

    def __init__(self, fallback_time_setting: TimeSettingType) -> None:
        """
        Initialize with a fallback time setting.

        Parameters
        ----------
        fallback_time_setting : TimeSettingType
            Fallback time setting to use if cached data is unavailable.
        """
        self.fallback_time_setting = resolve_time_setting(fallback_time_setting)

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Retrieve the timebase from cached data, or use the fallback if unavailable.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        if params.cache_data is not None:
            # set timebase to be the timebase of cached data
            try:
                times = params.cache_data["time"].to_numpy()
                # Check if the timebase is in ms instead of s
                if times[-1] > config(params.tokamak).MAX_SHOT_TIME:
                    times /= 1000  # [ms] -> [s]
                return times
            except KeyError:
                params.logger.warning(
                    "[Shot %s]: Shot constructor was passed data but no timebase.",
                    params.shot_id,
                )
                params.logger.debug(
                    "[Shot %s]: %s", params.shot_id, traceback.format_exc()
                )
        return self.fallback_time_setting.get_times(params)


class EfitTimeSetting(TimeSetting):
    """
    Time setting for using the EFIT timebase.
    """

    def __init__(self):
        """
        Initialize with tokamak-specific overrides for EFIT.
        """
        self.tokamak_overrides = {Tokamak.CMOD: self.cmod_times}

    def cmod_times(self, params: TimeSettingParams):
        """
        Retrieve the EFIT timebase for the CMOD tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        efit_tree_name = params.mds_conn.get_tree_name_of_nickname("_efit_tree")
        if efit_tree_name == "analysis":
            try:
                return params.mds_conn.get_data(
                    r"\analysis::efit_aeqdsk:time",
                    tree_name="_efit_tree",
                    astype="float64",
                )
            except mdsExceptions.MdsException:
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
        """
        Abstract method for retrieving EFIT timebase.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
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
        """
        Initialize with minimum current and duration values.

        Parameters
        ----------
        minimum_ip : float, optional
            Minimum current in amps (default is 400,000 A).
        minimum_duration : float, optional
            Minimum duration in seconds (default is 0.1 seconds).
        """
        self.tokamak_overrides = {Tokamak.D3D: self.d3d_times}
        self.minimum_ip = minimum_ip
        self.minimum_duration = minimum_duration

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Abstract method for retrieving disruption timebase.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        raise ValueError("Disruption time setting not implemented")

    def d3d_times(self, params: TimeSettingParams):
        """
        Retrieve disruption timebase for the D3D tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
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
        return times

    @classmethod
    def _get_end_of_shot(cls, signal, signal_time, threshold=1.0e5):
        """
        Calculate the end of shot based on signal and threshold.

        Parameters
        ----------
        signal : np.ndarray
            Signal array representing the current.
        signal_time : np.ndarray
            Time array corresponding to the signal.
        threshold : float
            Threshold value for determining shot end.

        Returns
        -------
        duration : float
            Duration of the shot.
        signal_max : float
            Maximum signal value.
        """
        duration = 0
        signal_max = 0
        if threshold < 0:
            raise Warning("Threshold is negative.")
        base_indices = np.where(signal_time <= 0.0)
        baseline = np.mean(signal[base_indices]) if len(base_indices) > 0 else 0
        signal -= baseline
        # Check if there was a finite signal otherwise consider the shot a "no plasma" shot
        finite_indices = np.where((signal_time >= 0.0) & (np.abs(signal) > threshold))
        if len(finite_indices) == 0:
            return duration, signal_max
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
    """
    Time setting for using the timebase of the plasma current.
    """

    def __init__(self):
        """
        Initialize with tokamak-specific overrides.
        """
        self.tokamak_overrides = {Tokamak.D3D: self.d3d_times}

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Abstract method for retrieving the Ip timebase.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        raise ValueError("Ip time setting not implemented")

    def d3d_times(self, params: TimeSettingParams):
        """
        Retrieve the Ip timebase for the D3D tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        (ip_time,) = params.mds_conn.get_dims(
            f"ptdata('ip', {params.shot_id})", tree_name="d3d"
        )
        return ip_time


class SignalTimeSetting(TimeSetting):
    """
    Time setting for using the timebase of a specific signal.
    """

    def __init__(self, tree_name: str, signal_path: str):
        """
        Initialize with the tree name and signal path.

        Parameters
        ----------
        tree_name : str
            Name of the tree containing the signal.
        signal_path : str
            Path to the signal within the tree.
        """
        self.tree_name = tree_name
        self.signal_path = signal_path

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Retrieve the timebase for the specified signal.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        try:
            (signal_time,) = params.mds_conn.get_dims(
                self.signal_path, tree_name=self.tree_name, astype="float64"
            )
            return signal_time
        except mdsExceptions.MdsException:
            params.logger.error(
                "Failed to set up timebase for signal %s", self.signal_path
            )
            raise


# --8<-- [start:time_setting_dict]
_time_setting_mappings: Dict[str, TimeSetting] = {
    "efit": EfitTimeSetting(),
    "disruption": DisruptionTimeSetting(),
    "disruption_warning": {
        Tokamak.CMOD: EfitTimeSetting(),
        Tokamak.D3D: DisruptionTimeSetting(),
    },
    "ip": IpTimeSetting(),
}
# --8<-- [end:time_setting_dict]


def resolve_time_setting(
    time_setting: TimeSettingType,
) -> TimeSetting:
    """
    Resolve a time setting to a TimeSetting object.

    Parameters
    ----------
    time_setting : TimeSettingType
        The time setting, which can be a string, list, or TimeSetting instance.

    Returns
    -------
    TimeSetting
        The resolved TimeSetting object.
    """
    if isinstance(time_setting, TimeSetting):
        return time_setting

    if isinstance(time_setting, str):
        time_setting_object = _time_setting_mappings.get(time_setting, None)
        if time_setting_object is not None:
            return time_setting_object

    if isinstance(time_setting, (list, np.ndarray)):
        return ListTimeSetting(time_setting)

    if isinstance(time_setting, pd.Series):
        return ListTimeSetting(time_setting.to_numpy())

    if isinstance(time_setting, dict):
        return TimeSettingDict(time_setting)

    raise ValueError("Invalid time setting")
