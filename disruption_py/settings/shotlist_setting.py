#!/usr/bin/env python3

from abc import ABC, abstractmethod
from dataclasses import dataclass
from importlib import resources
from logging import Logger
from typing import Dict, List, Type, Union

import numpy as np
import pandas as pd

import disruption_py.data
from disruption_py.io.sql import ShotDatabase
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.machine.tokamak import Tokamak


@dataclass
class ShotlistSettingParams:
    """Params passed by disruption_py to _get_shotlist() method.

    Attributes
    ----------
    database : ShotDatabase
        Database object to use for getting shotlist.
        A different database connection is used by each process.
        Defaults to logbook.
    tokamak : Tokamak
        The tokamak that is data is being retrieved for.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    database: ShotDatabase
    tokamak: Tokamak
    logger: Logger


class ShotlistSetting(ABC):
    """ShotlistSetting abstract class that should be inherited by all shotlist setting classes."""

    def get_shotlist(self, params: ShotlistSettingParams) -> List:
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_shotlist(params)

    @abstractmethod
    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        """Abstract method implemented by subclasses to get shotlist for the given setting params.

        Parameters
        ----------
        params : ShotlistSettingParams
            Params that can be used to determine shotlist.
        """
        pass


class IncludedShotlistSetting(ShotlistSetting):
    """Use the shotlist from one of the provided data files.

    Directly passing a key from the _get_shotlist_setting_mappings dictionary as a string will
    automatically create a new IncludedShotlistSetting object with that data_file_name.

    Parameters
    ----------
    data_file_name : str
        The name of the datafile that should be used to retrieve the shotlist.
    """

    def __init__(self, data_file_name: str) -> List:
        with resources.path(disruption_py.data, data_file_name) as data_file:
            df = pd.read_csv(data_file, header=None)
            lst = df.values[:, 0].tolist()
            self.shotlist = lst

    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        return self.shotlist


class FileShotlistSetting(ShotlistSetting):
    """Use a shotlist from the provided file path, this may be any file readable by pandas read_csv.

    Directly passing a file path as a string to the shotlist setting with the file name suffixed by txt or csv
    will automatically create a new FileShotlistSetting object with that file path.

    Parameters
    ----------
    file_path : str
        The file path of the file that should be used for retrieving the shotlist.
    column_index : int
        The index of the column that should be read. For text files, this should be 0. Defaults to 0.
    """

    def __init__(self, file_path, column_index=0):
        self.shotlist = (
            pd.read_csv(file_path, header=None).iloc[:, column_index].values.tolist()
        )

    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        return self.shotlist


class DatabaseShotlistSetting(ShotlistSetting):
    """Use an sql query of the database to retrieve the shotlist.

    Parameters
    ----------
    sql_query : str
        The sql query that should be used for retrieving shotlist.
    use_pandas : bool
        Whether Pandas should be used to do the sql query. Defaults to true.
    """

    def __init__(self, sql_query, use_pandas=True):
        self.sql_query = sql_query
        self.use_pandas = use_pandas

    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        if self.use_pandas:
            query_result_df = params.database.query(
                query=self.sql_query, use_pandas=True
            )
            return query_result_df.iloc[:, 0].tolist()
        else:
            query_result = params.database.query(query=self.sql_query, use_pandas=False)
            return [row[0] for row in query_result]


# --8<-- [start:get_shotlist_setting_dict]
_get_shotlist_setting_mappings: Dict[str, ShotlistSetting] = {
    "d3d_paper_shotlist": IncludedShotlistSetting("paper_shotlist.txt"),
    "d3d_train_disr": IncludedShotlistSetting("train_disr.txt"),
    "d3d_train_nondisr": IncludedShotlistSetting("train_nondisr.txt"),
    "cmod_test": IncludedShotlistSetting("cmod_test.txt"),
    "cmod_non_disruptions_ids_not_blacklist": IncludedShotlistSetting(
        "cmod_non_disruptions_ids_not_blacklist.txt"
    ),
    "cmod_non_disruptions_ids_not_blacklist_mini": IncludedShotlistSetting(
        "cmod_non_disruptions_ids_not_blacklist_mini.txt"
    ),
}
# --8<-- [end:get_shotlist_setting_dict]

# --8<-- [start:file_suffix_to_shotlist_setting_dict]
_file_suffix_to_shotlist_setting: Dict[str, Type[ShotlistSetting]] = {
    ".txt": FileShotlistSetting,
    ".csv": FileShotlistSetting,
}
# --8<-- [end:file_suffix_to_shotlist_setting_dict]

ShotlistSettingType = Union[
    "ShotlistSetting",
    int,
    str,
    Dict[Tokamak, "ShotlistSettingType"],
    List["ShotlistSettingType"],
]


def shotlist_setting_runner(shotlist_setting, params: ShotlistSettingParams):
    """
    Retrieve list of shot ids for the given shotlist setting.
    """
    if isinstance(shotlist_setting, ShotlistSetting):
        return shotlist_setting.get_shotlist(params)

    if isinstance(shotlist_setting, int) or (
        isinstance(shotlist_setting, str) and shotlist_setting.isdigit()
    ):
        return [shotlist_setting]

    if isinstance(shotlist_setting, np.ndarray):
        return shotlist_setting

    if isinstance(shotlist_setting, str):
        shotlist_setting_object = _get_shotlist_setting_mappings.get(
            shotlist_setting, None
        )
        if shotlist_setting_object is not None:
            return shotlist_setting_object.get_shotlist(params)

    if isinstance(shotlist_setting, str):
        # assume that it is a file path
        for suffix, shotlist_setting_type in _file_suffix_to_shotlist_setting.items():
            if shotlist_setting.endswith(suffix):
                return shotlist_setting_type(shotlist_setting).get_shotlist(params)

    if isinstance(shotlist_setting, dict):
        shotlist_setting = {
            map_string_to_enum(tokamak, Tokamak): shotlist_setting_mapping
            for tokamak, shotlist_setting_mapping in shotlist_setting.items()
        }
        chosen_setting = shotlist_setting.get(params.tokamak, None)
        if chosen_setting is not None:
            return shotlist_setting_runner(chosen_setting, params)

    if isinstance(shotlist_setting, list):
        all_results = []
        for setting in shotlist_setting:
            sub_result = shotlist_setting_runner(setting, params)
            if sub_result is not None:
                all_results.append(sub_result)

        return [shot_id for sub_list in all_results for shot_id in sub_list]

    raise ValueError("Invalid shot id setting")
