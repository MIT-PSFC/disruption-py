#!/usr/bin/env python3

"""
This module defines classes for time settings, used to manage the timebase for
retrieving data in disruption_py for various tokamaks and shot configurations.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Union

import numpy as np
from loguru import logger
from MDSplus import mdsExceptions

from disruption_py.config import config
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.core.utils.misc import shot_msg_patch
from disruption_py.inout.mds import MDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.east.util import EastUtilMethods
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
    database : ShotDatabase
        Database object with connection to the SQL database.
    disruption_time : float
        Time when the shot disrupted in seconds (or None if no disruption occurred).
    tokamak : Tokamak
        Tokamak for which the time setting is applied.
    """

    shot_id: int
    mds_conn: MDSConnection
    database: ShotDatabase
    disruption_time: float
    tokamak: Tokamak

    def __post_init__(self):
        self.logger = shot_msg_patch(logger, self.shot_id)

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


TimeSettingType = Union[
    "TimeSetting", str, Dict[Tokamak, "TimeSettingType"], List["TimeSettingType"]
]


def _postprocess(times: np.ndarray, units: str = "") -> np.ndarray:
    """
    Return a deduplicated, sorted, rescaled time array.

    Parameters
    ----------
    times : np.ndarray
        Array of times.
    units : str
        Unit of measure for the array of times.

    Returns
    -------
    np.ndarray
        Post-processed array of times.
    """
    times = np.unique(times)
    units = units.lower().strip()
    if units == "ms":
        return times * 1e-3
    if units == "us":
        return times * 1e-6
    return times


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
                return _postprocess(
                    times=self.tokamak_overrides[params.tokamak](params)
                )
        return _postprocess(times=self._get_times(params))

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
        params.logger.warning(
            "No time setting for tokamak {tokamak}", tokamak=params.tokamak
        )
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
        times : list, np.ndarray
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


class EfitTimeSetting(TimeSetting):
    """
    Time setting for using the EFIT timebase.
    """

    def _get_times(self, params: TimeSettingParams):
        """
        Retrieve the EFIT timebase for the tested tokamaks.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        (efit_time,) = params.mds_conn.get_dims(
            r"\efit_aeqdsk:ali", tree_name="_efit_tree", astype="float64"
        )
        efit_time_unit = params.mds_conn.get_data(
            r"units_of(dim_of(\efit_aeqdsk:ali))", tree_name="_efit_tree", astype="str"
        )
        if efit_time_unit != "s":
            params.logger.verbose(
                "Failed to get the time units of EFIT tree '{tree}', assuming seconds.",
                tree=params.mds_conn.get_tree_name_of_nickname("_efit_tree"),
            )
        return _postprocess(times=efit_time, units=efit_time_unit)


class DisruptionTimeSetting(TimeSetting):
    """
    Time setting for using the disruption timebase.

    The disruption timebase is the timebase of the shot that was disrupted.
    """

    # Disruption Variables
    DT_BEFORE_DISRUPTION_D3D = 0.002
    DURATION_BEFORE_DISRUPTION_D3D = 0.10
    DT_BEFORE_DISRUPTION_EAST = 0.010
    DURATION_BEFORE_DISRUPTION_EAST = 0.25

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
        self.tokamak_overrides = {
            Tokamak.D3D: self.d3d_times,
            Tokamak.EAST: self.east_times,
        }
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
        baseline = np.mean(raw_ip[:10])
        ip = raw_ip - baseline
        duration, ip_max = self._get_end_of_shot(ip, ip_time, 100e3)
        if duration < self.minimum_duration or np.abs(ip_max) < self.minimum_ip:
            raise NotImplementedError()

        times = np.arange(0.100, duration + config(params.tokamak).time_const, 0.025)
        if params.disrupted:
            additional_times = np.arange(
                params.disruption_time - self.DURATION_BEFORE_DISRUPTION_D3D,
                params.disruption_time + config(params.tokamak).time_const,
                self.DT_BEFORE_DISRUPTION_D3D,
            )
            times = times[
                np.where(
                    times
                    < (
                        params.disruption_time
                        - self.DURATION_BEFORE_DISRUPTION_D3D
                        - config(params.tokamak).time_const
                    )
                )
            ]
            times = np.concatenate((times, additional_times))
        return times

    def east_times(self, params: TimeSettingParams):
        """
        Retrieve disruption timebase for the EAST tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        ip, ip_time = EastUtilMethods.retrieve_ip(params.mds_conn, params.shot_id)
        # For EAST, minimum_ip = 200e3 [A], minimum_duration = 0.6 [s]
        duration, ip_max = self._get_end_of_shot(ip, ip_time, self.minimum_ip)
        if duration < self.minimum_duration or np.abs(ip_max) < self.minimum_ip:
            raise NotImplementedError()

        times = np.arange(0.200, duration + config(params.tokamak).time_const, 0.1)
        if params.disrupted:
            additional_times = np.arange(
                params.disruption_time - self.DURATION_BEFORE_DISRUPTION_EAST,
                params.disruption_time + config(params.tokamak).time_const,
                self.DT_BEFORE_DISRUPTION_EAST,
            )
            times = times[
                np.where(
                    times
                    < (
                        params.disruption_time
                        - self.DURATION_BEFORE_DISRUPTION_EAST
                        - config(params.tokamak).time_const
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
        (base_indices,) = np.where(signal_time <= 0.0)
        baseline = np.mean(signal[base_indices]) if len(base_indices) > 0 else 0
        signal -= baseline
        # Check if there was a finite signal otherwise consider the shot a "no plasma" shot
        (finite_indices,) = np.where(
            (signal_time >= 0.0) & (np.abs(signal) > threshold)
        )
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
        (valid_indices,) = np.where(
            (polarized_signal >= threshold) & (signal_time > 0.0)
        )
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
        self.tokamak_overrides = {
            Tokamak.CMOD: self.cmod_times,
            Tokamak.D3D: self.d3d_times,
            Tokamak.EAST: self.east_times,
        }

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

    def cmod_times(self, params: TimeSettingParams):
        """
        Retrieve the Ip timebase for the CMOD tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        (ip_time,) = params.mds_conn.get_dims(r"\ip", tree_name="magnetics")
        return ip_time

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
            f"ptdata('ip', {params.shot_id})", tree_name=None
        )
        ip_time /= 1e3  # [ms] -> [s]
        return ip_time

    def east_times(self, params: TimeSettingParams):
        """
        Retrieve the Ip timebase for the EAST tokamak.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        (ip_time,) = params.mds_conn.get_dims(r"\pcrl01", tree_name="pcs_east")
        # For shots before year 2014, the PCRL01 timebase needs to be shifted
        # by 17.0 ms
        if params.shot_id < 44432:
            ip_time -= 0.0170
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
        self.tokamak_overrides = {
            Tokamak.D3D: self._get_times_d3d_wrapper,
        }

    def _get_times_d3d_wrapper(self, params: TimeSettingParams) -> np.ndarray:
        """
        Wrapper function for the DIII-D tokamak
        - Allow using tree_name="ptdata" to call PTDATA signals
        - Convert time unit from ms to s
        """
        if self.tree_name.lower() == "ptdata":
            self.signal_path = f'ptdata("{self.signal_path}", {params.shot_id})'
            self.tree_name = None
        return self._get_times(params)

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
        except mdsExceptions.MdsException:
            params.logger.error(
                "Failed to set up timebase for signal {signal_path}",
                signal_path=self.signal_path,
            )
            raise
        signal_unit = params.mds_conn.get_data(
            f"units_of(dim_of({self.signal_path}))",
            tree_name=self.tree_name,
            astype="str",
        )
        if (
            not signal_unit.strip()
            and self.signal_path.startswith("ptdata")
            and params.tokamak == Tokamak.D3D
        ):
            # timebase of PTDATA signal defaults to [ms]
            signal_unit = "ms"
        elif signal_unit != "s":
            params.logger.warning(
                "Failed to get the time unit of signal `{path}`, assuming seconds.",
                path=self.signal_path,
            )
        return _postprocess(times=signal_time, units=signal_unit)


class SharedTimeSetting(TimeSetting):
    """
    Time setting for using a "shared" time among several time settings.
    The first time setting in the list is used to determine the timebase, while the other
    time settings in the list are used to select the starting and ending points, so that
    the final time base focuses on a time window shared among all the given items.
    """

    def __init__(self, time_setting_list: list[TimeSetting]):
        """
        Initialize with a list of time settings.

        Parameters
        ----------
        time_setting_list : list[TimeSetting]
            List of time settings.
        """
        if len(time_setting_list) < 2:
            raise ValueError("SharedTimeSetting requires at least 2 time settings.")
        self.time_setting_list = time_setting_list

    def _get_times(self, params: TimeSettingParams) -> np.ndarray:
        """
        Retrieve the shared timebase for the given settings.

        Parameters
        ----------
        params : TimeSettingParams
            Parameters needed to retrieve the timebase.

        Returns
        -------
        np.ndarray
            Array of times in the timebase.
        """
        times, *others = [ts.get_times(params) for ts in self.time_setting_list]
        tmin = np.max([np.min(t) for t in others])
        tmax = np.min([np.max(t) for t in others])
        return times[np.where((times >= tmin) & (times <= tmax))]


# --8<-- [start:time_setting_dict]
_time_setting_mappings: Dict[str, TimeSetting] = {
    "efit": EfitTimeSetting(),
    "disruption": {
        Tokamak.D3D: DisruptionTimeSetting(),
        Tokamak.EAST: DisruptionTimeSetting(minimum_ip=200e3, minimum_duration=0.6),
    },
    "disruption_warning": {
        Tokamak.CMOD: EfitTimeSetting(),
        Tokamak.D3D: DisruptionTimeSetting(),
        Tokamak.EAST: DisruptionTimeSetting(minimum_ip=200e3, minimum_duration=0.6),
    },
    "ip": IpTimeSetting(),
    "ip_efit": SharedTimeSetting([IpTimeSetting(), EfitTimeSetting()]),
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
        time_setting_object = _time_setting_mappings.get(time_setting)
        if time_setting_object is not None:
            return time_setting_object

    if isinstance(time_setting, list):
        if all(isinstance(ts, TimeSetting) for ts in time_setting):
            return SharedTimeSetting(time_setting)
        return ListTimeSetting(time_setting)

    if isinstance(time_setting, np.ndarray):
        return ListTimeSetting(time_setting)

    if isinstance(time_setting, dict):
        return TimeSettingDict(time_setting)

    raise ValueError("Invalid time setting")
