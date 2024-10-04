#!/usr/bin/env python3

"""
Handles retrieving shotlists from various sources including lists, files, and SQL 
databases.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from importlib import resources
from logging import Logger
from typing import Dict, List, Type, Union

import numpy as np
import pandas as pd

import disruption_py.data
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class ShotlistSettingParams:
    """
    Params passed by disruption_py to _get_shotlist() method.

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
        """
        Retrieve the shotlist based on the provided parameters.

        Parameters
        ----------
        params : ShotlistSettingParams
            The parameters containing the database, tokamak, and logger used
            to determine the shotlist.

        Returns
        -------
        List
            A list of shot IDs retrieved.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_shotlist(params)

    @abstractmethod
    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        """
        Abstract method implemented by subclasses to get shotlist for the given setting params.

        Parameters
        ----------
        params : ShotlistSettingParams
            Params that can be used to determine shotlist.
        """


class FileShotlistSetting(ShotlistSetting):
    """
    Use `pandas.read_csv` to read a file, then extract and use values from any column.

    Directly passing a file path as a string to the shotlist setting with the file name suffixed
    by txt or csv will automatically create a new FileShotlistSetting object with that file path.

    Parameters
    ----------
    file_path : str
        The file path of the file that should be used for retrieving the shotlist.
    column_index : int
        The index of the column that should be read. Defaults to 0.
    **kwargs : Dict
        Optional keyword arguments dictionary to be passed to `pandas.read_csv`.
    """

    def __init__(self, file_path: str, column_index: int = 0, **kwargs: Dict):
        self.file_path = file_path
        self.column_index = column_index
        self.kwargs = kwargs
        self.shotlist = []

    def _get_shotlist(self, params: ShotlistSettingParams) -> List:
        if not self.shotlist:
            self.kwargs.setdefault("header", "infer")
            df = pd.read_csv(self.file_path, **self.kwargs)
            arr = df.values[:, self.column_index]
            self.shotlist = arr.astype(int).tolist()
        return self.shotlist


class IncludedShotlistSetting(FileShotlistSetting):
    """
    Use the shotlist from one of the provided data files.

    Directly passing a key from the _get_shotlist_setting_mappings dictionary as a string will
    automatically create a new IncludedShotlistSetting object with that file_name.

    Parameters
    ----------
    file_name : str
        The name of the datafile that should be used to retrieve the shotlist.
    **kwargs : Dict
        Optional keyword arguments dictionary to be passed to `FileShotlistSetting`.
    """

    def __init__(self, file_name: str, **kwargs: Dict):
        data = resources.files(disruption_py.data)
        file = data.joinpath(file_name)
        with resources.as_file(file) as file_path:
            super().__init__(file_path, **kwargs)


class DatabaseShotlistSetting(ShotlistSetting):
    """
    Use an sql query of the database to retrieve the shotlist.

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
        query_result = params.database.query(query=self.sql_query, use_pandas=False)
        return [row[0] for row in query_result]


# --8<-- [start:get_shotlist_setting_dict]
_get_shotlist_setting_mappings: Dict[str, ShotlistSetting] = {
    "disruption_warning": DatabaseShotlistSetting(
        "select distinct shot from disruption_warning"
    ),
    "plasmas": DatabaseShotlistSetting(
        """
        if exists (select * from information_schema.tables where table_name = 'summary')
        begin
            select distinct shot from summary where ipmax > 100e3 and pulse_length > 0.1;
        end
        else if exists (select * from information_schema.tables where table_name = 'summaries')
        begin
            select distinct shot from summaries where ipmax > 100e3 and pulse_length > 0.1;
        end
        """
    ),
    "cmod_ufo": IncludedShotlistSetting("cmod_ufo.csv"),
    "cmod_vde": IncludedShotlistSetting("cmod_vde.csv"),
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
        # Do not immediately return the list because it may be multidimensional
        # and would need to be handled as such below
        shotlist_setting = shotlist_setting.get_shotlist(params)

    if isinstance(shotlist_setting, (int, np.int64)) or (
        isinstance(shotlist_setting, str) and shotlist_setting.isdigit()
    ):
        return [shotlist_setting]

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

    if isinstance(shotlist_setting, (list, np.ndarray)):
        all_results = []
        for setting in shotlist_setting:
            sub_result = shotlist_setting_runner(setting, params)
            if sub_result is not None:
                all_results.append(sub_result)

        return [shot_id for sub_list in all_results for shot_id in sub_list]

    raise ValueError("Invalid shot id setting")
