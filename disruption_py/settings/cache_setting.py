#!/usr/bin/env python3

"""
This module provides classes for managing and retrieving cached data.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, Type, Union

import xarray as xr
from loguru import logger

from disruption_py.core.utils.misc import shot_log_patch
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class CacheSettingParams:
    """
    Params passed by disruption_py to _get_cache_data() method.

    Attributes
    ----------
    shot_id : int
        Shot Id for which to get cached data.
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
        self.logger = shot_log_patch(logger, self.shot_id)


CacheSettingType = Union["CacheSetting", xr.Dataset, str]


class CacheSetting(ABC):
    """
    CacheSetting abstract class that should be inherited by all cache setting classes.

    Subclasses must implement the `_get_cache_data` method to define how cached data is retrieved.
    """

    def get_cache_data(self, params: CacheSettingParams) -> xr.Dataset:
        """
        Return cached data, using tokamak-specific overrides if available.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_cache_data(params)

    @abstractmethod
    def _get_cache_data(self, params: CacheSettingParams) -> xr.Dataset:
        """
        Abstract method implemented by subclasses to get cached data for a
        given set of params as an xarray dataset.

        Parameters
        ----------
        params : CacheSettingParams
            Params that can be used to determine and retrieve cached data.

        Returns
        -------
        xr.Dataset
            xarray dataset containing cached data.
        """


class SQLCacheSetting(CacheSetting):
    """
    Cache setting for retrieving data from the SQL database.
    """

    def _get_cache_data(self, params: CacheSettingParams) -> xr.Dataset:
        df = params.database.get_shots_data(shotlist=[params.shot_id])
        old = df.shape[0]
        params.logger.debug(
            "Retrieved cache from SQL: {old} rows",
            old=old,
        )
        df.drop_duplicates(subset=["shot", "time"], inplace=True)
        new = df.shape[0]
        if new < old:
            params.logger.debug(
                "Pruned cache: {diff} rows",
                diff=old - new,
            )
        ds = xr.Dataset.from_dataframe(df.set_index(["shot", "time"]))
        return ds


class DatasetCacheSetting(CacheSetting):
    """
    Cache setting for retrieving data from an xarray Dataset.

    Parameters
    ----------
    cache_setting : xr.Dataset or str
        The Dataset to use as the cached data, or its file path.
    """

    def __init__(self, cache_setting: xr.Dataset | str):
        if isinstance(cache_setting, str):
            cache_data = xr.open_dataset(cache_setting)
        elif isinstance(cache_setting, xr.Dataset):
            cache_data = cache_setting
        else:
            raise ValueError(f"Invalid cache setting: {cache_setting}")
        self.cache_data = cache_data

    def _get_cache_data(self, params: CacheSettingParams) -> xr.Dataset:
        return self.cache_data


# --8<-- [start:cache_setting_dict]
_cache_setting_mappings: Dict[str, CacheSetting] = {
    "sql": SQLCacheSetting,
}
# --8<-- [end:cache_setting_dict]

_file_suffix_to_cache_setting: Dict[str, Type[CacheSetting]] = {
    ".cdf": DatasetCacheSetting,
    ".h5": DatasetCacheSetting,
    ".hdf5": DatasetCacheSetting,
    ".nc": DatasetCacheSetting,
}


def resolve_cache_setting(
    cache_setting: CacheSettingType,
) -> CacheSetting | None:
    """
    Resolve the cache setting to a CacheSetting instance.

    Parameters
    ----------
    cache_setting : CacheSettingType
        The cache setting to resolve. This can be an instance of CacheSetting,
        a string representing a cache setting type, or an xarray Dataset.

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
        # shortcuts
        cache_setting_cls = _cache_setting_mappings.get(cache_setting)
        if cache_setting_cls:
            return cache_setting_cls()
        # extensions
        for ext, cache_setting_cls in _file_suffix_to_cache_setting.items():
            if cache_setting.lower().endswith(ext):
                return cache_setting_cls(cache_setting)

    if isinstance(cache_setting, xr.Dataset):
        return DatasetCacheSetting(cache_setting)

    raise ValueError("Invalid cache setting")
