#!/usr/bin/env python3

from dataclasses import dataclass
from logging import Logger
from abc import ABC, abstractmethod
from typing import Dict, Union

from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.io.mds import MDSConnection
from disruption_py.io.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak

NicknameSettingType = Union[
    "NicknameSettingType", str, Dict[Tokamak, "NicknameSettingType"]
]


@dataclass
class NicknameSettingParams:
    """Params passed by disruption_py to nickname trees.

    Attributes
    ----------
    shot_id : int
        the shot id for which to resolve nicknames.
    mds_conn: MDSConnection
        MDSConnection object to access MDSPlus data.
    database : ShotDatabase
        Database connection object for tokamak sql database.
        A different database connection is used by each thread/process.
    disruption_time : float
        The time of the disruption in seconds.
    tokamak : Tokamak
        The tokamak for which results are being output.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    shot_id: int
    mds_conn: MDSConnection
    database: ShotDatabase
    disruption_time: float
    tokamak: Tokamak
    logger: Logger


class NicknameSetting(ABC):
    """
    A setting for getting tree nicknames.
    """

    def resolve_nickname_func(self, params: NicknameSettingParams):
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return lambda: self.tokamak_overrides[params.tokamak](params)
        return lambda: self._resolve_nickname(params)

    @abstractmethod
    def _resolve_nickname(self, params: NicknameSettingParams) -> str:
        """Abstract method implemented by subclasses to determine an MDSplus tree name.

        Parameters
        ----------
        params : NicknameSettingParams
            Params that can be used to determine and retrieve the MDSplus tree name.

        Returns
        -------
        str
            The MDSplus tree name.
        """


class NicknameSettingDict(NicknameSetting):
    """
    Utility class that is automatically used when a dictionary is passed the nickname setting.

    Parameters
    ----------
    nickname_setting_dict : dict[Tokamak, NicknameSettingType]
        A dictionary mapping tokamaks to a nickname setting, e.g. `{'cmod': 'efit'}`.
        Any other option passable as a value in the `nickname_setting_dict` dictionary parameter
        in `RetrievalSettings` may be used.
    """

    def __init__(self, nickname_setting_dict: Dict[Tokamak, NicknameSettingType]):
        resolved_nickname_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_nickname_setting(
                individual_setting
            )
            for tokamak, individual_setting in nickname_setting_dict.items()
        }
        self.resolved_nickname_setting_dict = resolved_nickname_setting_dict

    def _resolve_nickname(self, params: NicknameSettingParams) -> str:
        chosen_setting = self.resolved_nickname_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_times(params)
        raise NotImplementedError(
            f"{self.__class__.__name__} is not implemented for tokamak {params.tokamak}."
        )


class StaticNicknameSetting(NicknameSetting):
    def __init__(self, tree_name: str):
        self.tree_name = tree_name

    def _resolve_nickname(self, params: NicknameSettingParams) -> str:
        return self.tree_name


class DefaultNicknameSetting(NicknameSetting):
    """
    A setting to resolve the '_efit_tree' nickname to the default EFIT tree.
    """

    def __init__(self):
        self.tokamak_overrides = {
            Tokamak.CMOD: StaticNicknameSetting("analysis")._resolve_nickname,
            Tokamak.D3D: StaticNicknameSetting("efit01")._resolve_nickname,
        }

    def _resolve_nickname(self, params: NicknameSettingParams) -> str:
        raise NotImplementedError(
            f"{self.__class__.__name__} is not implemented for tokamak {params.tokamak}."
        )


class DisruptionNicknameSetting(NicknameSetting):
    """
    A setting to resolve the '_efit_tree' nickname to the disruption EFIT tree.
    """

    def __init__(self):
        self.tokamak_overrides = {
            Tokamak.CMOD: self._cmod_nickname,
            Tokamak.D3D: self._d3d_nickname,
        }

    def _d3d_nickname(self, params: NicknameSettingParams) -> str:
        if params.disruption_time is None:
            return DefaultNicknameSetting(params)._resolve_nickname()
        efit_trees = params.database.query(
            "select tree from code_rundb.dbo.plasmas where "
            f"shot = {params.shot_id} and runtag = 'DIS' and deleted = 0 order by idx",
            use_pandas=False,
        )
        if len(efit_trees) == 0:
            return DefaultNicknameSetting(params)._resolve_nickname()
        efit_tree = efit_trees[-1][0]
        return efit_tree

    def _cmod_nickname(self, params: NicknameSettingParams) -> str:
        if params.disruption_time is None:
            return DefaultNicknameSetting(params)._resolve_nickname()
        return "efit18"

    def _resolve_nickname(self, params: NicknameSettingParams) -> str:
        raise NotImplementedError(
            f"{self.__class__.__name__} is not implemented for tokamak {params.tokamak}."
        )


# --8<-- [start:nickname_setting_keys]
_nickname_setting_mappings: Dict[str, NicknameSetting] = {
    "default": DefaultNicknameSetting(),
    "disruption": DisruptionNicknameSetting(),
    # deprecated
    "analysis": DefaultNicknameSetting(),
    "disruption_warning": DisruptionNicknameSetting(),
}
# --8<-- [end:nickname_setting_keys]


def resolve_nickname_setting(nickname_setting: NicknameSettingType) -> NicknameSetting:
    """
    Resolve a nickname setting to a nickname function.
    """
    if isinstance(nickname_setting, NicknameSetting):
        return nickname_setting
    if isinstance(nickname_setting, dict):
        return NicknameSettingDict(nickname_setting)
    if (
        isinstance(nickname_setting, str)
        and nickname_setting in _nickname_setting_mappings
    ):
        return _nickname_setting_mappings[nickname_setting]

    if isinstance(nickname_setting, str):
        return StaticNicknameSetting(nickname_setting)

    raise ValueError(f"Invalid nickname setting type {type(nickname_setting)}.")
