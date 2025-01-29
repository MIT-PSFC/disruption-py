#!/usr/bin/env python3

"""
This module provides classes for managing and retrieving cached data from various
sources, including SQL databases and Pandas DataFrames.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, Type, Union

import pandas as pd
from loguru import logger

from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.core.utils.misc import shot_log_msg
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class CacheSettingParams:
    """
    Params passed by disruption_py to _get_cache_data() method.

    Attributes
    ----------
    shot_id : int
        Shot Id for which to get cached data. Defaults to logbook.
    database : ShotDatabase
        Database object to use for getting cached data.
        A different database connection is used by each thread/process.
    tokamak : Tokamak
        The tokamak being run.
    """

    shot_id: int
    database: ShotDatabase
    tokamak: Tokamak

    def __post_init__(self):
        self.logger = logger.patch(
            lambda record: record.update(
                message=shot_log_msg(self.shot_id, record["message"])
            )
        )


CacheSettingType = Union[
    "CacheSetting", str, pd.DataFrame, Dict[Tokamak, "CacheSettingType"]
]


class CacheSetting(ABC):
    """
    CacheSetting abstract class that should be inherited by all cache setting classes.

    Subclasses must implement the `_get_cache_data` method to define how cached data is retrieved.
    """

    def get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        """
        Return cached data, using tokamak-specific overrides if available.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_cache_data(params)

    @abstractmethod
    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        """
        Abstract method implemented by subclasses to get cached data for a
        given set of params as a Pandas dataframe.

        Parameters
        ----------
        params : CacheSettingParams
            Params that can be used to determine and retrieve cached data.

        Returns
        -------
        pd.DataFrame
            Pandas dataframe containing cached data.
        """


class CacheSettingDict(CacheSetting):
    """
    A class that resolves and manages cache settings for Tokamak instances.

    This class takes a dictionary of cache settings and resolves them for easy access.

    Attributes
    ----------
    resolved_cache_setting_dict : Dict[Tokamak, CacheSettingType]
        A mapping of Tokamak instances to their resolved cache settings.

    Parameters
    ----------
    cache_setting_dict : Dict[Tokamak, CacheSettingType]
        A dictionary of initial cache settings for each Tokamak.

    Methods
    -------
    _get_cache_data(params: CacheSettingParams) -> pd.DataFrame
        Retrieves cache data for the specified Tokamak.
    """

    def __init__(self, cache_setting_dict: Dict[Tokamak, CacheSettingType]):
        resolved_cache_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_cache_setting(
                individual_setting
            )
            for tokamak, individual_setting in cache_setting_dict.items()
        }
        self.resolved_cache_setting_dict = resolved_cache_setting_dict

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        chosen_setting = self.resolved_cache_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_cache_data(params)
        params.logger.warning(
            "No cache setting for tokamak {tokamak}", tokamak=params.tokamak
        )
        return None


class SQLCacheSetting(CacheSetting):
    """Cache setting for retrieving data from SQL database."""

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        params.logger.info("retrieving sql data")
        return params.database.get_shots_data(shotlist=[params.shot_id])


class DFCacheSetting(CacheSetting):
    """
    Cache setting for retrieving data from a Pandas DataFrame.

    Parameters
    ----------
    cache_data : pd.DataFrame
        The DataFrame to use as the cached data.
    """

    def __init__(self, cache_data: pd.DataFrame):
        self.cache_data = cache_data

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        return self.cache_data


class HDF5CacheSetting(CacheSetting):
    """
    Cache setting for retrieving data from an HDF5 file.

    Parameters
    ----------
    cache_file : str
        The path to the HDF5 file to use as the cached data.
    """

    def __init__(self, cache_file: str):
        self.cache_file = cache_file

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        # Some shots may not have been cached, so the key will not exist, and
        # it throws an error. If there is no cached data, return None.
        try:
            return pd.read_hdf(self.cache_file, key=f"df_{params.shot_id}")
        except KeyError:
            return None


class CSVCacheSetting(CacheSetting):
    """
    Cache setting for retrieving data from a CSV file.

    Parameters
    ----------
    cache_file : str
        The path to the CSV file to use as the cached data.
    """

    def __init__(self, cache_file: str):
        self.cache_file = cache_file

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        return pd.read_csv(self.cache_file)


# --8<-- [start:cache_setting_dict]
_cache_setting_mappings: Dict[str, CacheSetting] = {
    "sql": SQLCacheSetting(),
}
# --8<-- [end:cache_setting_dict]

_file_suffix_to_cache_setting: Dict[str, Type[CacheSetting]] = {
    ".h5": HDF5CacheSetting,
    ".hdf5": HDF5CacheSetting,
    ".csv": CSVCacheSetting,
}


def resolve_cache_setting(
    cache_setting: CacheSettingType,
) -> CacheSetting:
    """
    Resolve the cache setting to a CacheSetting instance.

    Parameters
    ----------
    cache_setting : CacheSettingType
        The cache setting to resolve. This can be an instance of CacheSetting,
        a string representing a cache setting type, a Pandas DataFrame, or a
        dictionary of cache settings.

    Returns
    -------
    CacheSetting
        An instance of CacheSetting corresponding to the provided cache setting.
    """
    if cache_setting is None:
        return None

    if isinstance(cache_setting, CacheSetting):
        return cache_setting

    if isinstance(cache_setting, str):
        cache_setting_object = _cache_setting_mappings.get(cache_setting, None)
        if cache_setting_object is not None:
            return cache_setting_object

        # assume that it is a file path
        for (
            suffix,
            cache_setting_type,
        ) in _file_suffix_to_cache_setting.items():
            if cache_setting.endswith(suffix):
                return cache_setting_type(cache_setting)

    if isinstance(cache_setting, pd.DataFrame):
        return DFCacheSetting(cache_setting)

    if isinstance(cache_setting, dict):
        return CacheSettingDict(cache_setting)

    raise ValueError("Invalid cache setting")
