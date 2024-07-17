#!/usr/bin/env python3

from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Dict, Union

import pandas as pd

from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.io.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class CacheSettingParams:
    """Params passed by disruption_py to _get_cache_data() method.

    Attributes
    ----------
    shot_id : int
        Shot Id for which to get cached data. Defaults to logbook.
    database : ShotDatabase
        Database object to use for getting cached data.
        A different database connection is used by each thread/process.
    tokamak : Tokamak
        The tokamak being run.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    shot_id: int
    database: ShotDatabase
    tokamak: Tokamak
    logger: Logger


CacheSettingType = Union[
    "CacheSetting", str, pd.DataFrame, Dict[Tokamak, "CacheSettingType"]
]


class CacheSetting(ABC):
    """CacheSetting abstract class that should be inherited by all cache setting classes."""

    def get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_cache_data(params)

    @abstractmethod
    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        """Abstract method implemented by subclasses to get cached data for a given set of params as a pandas dataframe.

        Parameters
        ----------
        params : CacheSettingParams
            Params that can be used to determine and retrieve cached data.

        Returns
        -------
        pd.DataFrame
            Pandas dataframe containing cached data.
        """
        pass


class CacheSettingDict(CacheSetting):
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
        else:
            params.logger.warning(f"No cache setting for tokamak {params.tokamak}")
            return None


class SQLCacheSetting(CacheSetting):
    """Cache setting for retrieving data from SQL database."""

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        params.logger.info(f"retrieving sql data for {params.shot_id}")
        return params.database.get_shots_data(shotlist=[params.shot_id])


class DFCacheSetting(CacheSetting):
    """Cache setting for retrieving data from a pandas dataframe.

    Parameters
    ----------
    cache_data : pd.DataFrame
        The dataframe to use as the cached data.
    """

    def __init__(self, cache_data: pd.DataFrame):
        self.cache_data = cache_data

    def _get_cache_data(self, params: CacheSettingParams) -> pd.DataFrame:
        return self.cache_data


# --8<-- [start:cache_setting_dict]
_cache_setting_mappings: Dict[str, CacheSetting] = {
    "sql": SQLCacheSetting(),
}
# --8<-- [end:cache_setting_dict]


def resolve_cache_setting(
    cache_setting: CacheSettingType,
) -> CacheSetting:
    if cache_setting is None:
        return None

    if isinstance(cache_setting, CacheSetting):
        return cache_setting

    if isinstance(cache_setting, str):
        cache_setting = _cache_setting_mappings.get(cache_setting, None)
        if cache_setting is not None:
            return cache_setting

    if isinstance(cache_setting, pd.DataFrame):
        return DFCacheSetting(cache_setting)

    if isinstance(cache_setting, dict):
        return CacheSettingDict(cache_setting)

    raise ValueError("Invalid cache setting")
