#!/usr/bin/env python3

from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Dict, Union

import pandas as pd

from disruption_py.io.sql import ShotDatabase
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.machine.tokamak import Tokamak


@dataclass
class InputSettingParams:
    """Params passed by disruption_py to _get_input_data() method.

    Attributes
    ----------
    shot_id : int
        Shot Id for which to get input data. Defaults to logbook.
    database : ShotDatabase
        Database object to use for getting input data.
        A different database connection is used by each thread/process.
    tokamak : Tokemak
        The tokamak being run.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    shot_id: int
    database: ShotDatabase
    tokamak: Tokamak
    logger: Logger


InputSettingType = Union[
    "InputSetting", str, pd.DataFrame, Dict[Tokamak, "InputSettingType"]
]


class InputSetting(ABC):
    """InputSetting abstract class that should be inherited by all input setting classes."""

    def get_input_data(self, params: InputSettingParams) -> pd.DataFrame:
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_input_data(params)

    @abstractmethod
    def _get_input_data(self, params: InputSettingParams) -> pd.DataFrame:
        """Abstract method implemented by subclasses to get input data for a given set of params as a pandas dataframe.

        Parameters
        ----------
        params : InputSettingParams
            Params that can be used to determine and retrieve input data.

        Returns
        -------
        pd.DataFrame
            Pandas dataframe containing input data.
        """
        pass


class InputSettingDict(InputSetting):
    def __init__(self, input_setting_dict: Dict[Tokamak, InputSettingType]):
        resolved_input_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_input_setting(
                individual_setting
            )
            for tokamak, individual_setting in input_setting_dict.items()
        }
        self.resolved_input_setting_dict = resolved_input_setting_dict

    def _get_input_data(self, params: InputSettingParams) -> pd.DataFrame:
        chosen_setting = self.resolved_input_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_input_data(params)
        else:
            params.logger.warning(f"No input setting for tokamak {params.tokamak}")
            return None


class SQLInputSetting(InputSetting):
    """Input setting for retrieving data from SQL database."""

    def _get_input_data(self, params: InputSettingParams) -> pd.DataFrame:
        params.logger.info(f"retrieving sql data for {params.shot_id}")
        return params.database.get_shots_data(shotlist=[params.shot_id])


class DFInputSetting(InputSetting):
    """Input setting for retrieving data from a pandas dataframe.

    Parameters
    ----------
    input_data : pd.DataFrame
        The dataframe to use as the input data.
    """

    def __init__(self, input_data: pd.DataFrame):
        self.input_data = input_data

    def _get_input_data(self, params: InputSettingParams) -> pd.DataFrame:
        return self.input_data


# --8<-- [start:input_setting_dict]
_input_setting_mappings: Dict[str, InputSetting] = {
    "sql": SQLInputSetting(),
}
# --8<-- [end:input_setting_dict]


def resolve_input_setting(
    input_setting: InputSettingType,
) -> InputSetting:
    if input_setting is None:
        return None

    if isinstance(input_setting, InputSetting):
        return input_setting

    if isinstance(input_setting, str):
        input_setting = _input_setting_mappings.get(input_setting, None)
        if input_setting is not None:
            return input_setting

    if isinstance(input_setting, pd.DataFrame):
        return DFInputSetting(input_setting)

    if isinstance(input_setting, dict):
        return InputSettingDict(input_setting)

    raise ValueError("Invalid input setting")
