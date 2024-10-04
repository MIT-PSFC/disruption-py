#!/usr/bin/env python3

"""
This module provides classes for handling nickname settings, which are used to 
resolve MDSplus tree names for various tokamaks.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Dict, Union

from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.inout.mds import MDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak

NicknameSettingType = Union[
    "NicknameSettingType", str, Dict[Tokamak, "NicknameSettingType"]
]


@dataclass
class NicknameSettingParams:
    """
    Parameters passed to nickname trees for resolving tree nicknames.

    Attributes
    ----------
    shot_id : int
        The shot ID for which to resolve nicknames.
    mds_conn : MDSConnection
        MDSConnection object for accessing MDSPlus data.
    database : ShotDatabase
        Database connection for querying tokamak shot data.
    disruption_time : float
        The time of the disruption in seconds.
    tokamak : Tokamak
        The tokamak for which results are being processed.
    logger : Logger
        Logger for logging relevant messages.
    """

    shot_id: int
    mds_conn: MDSConnection
    database: ShotDatabase
    disruption_time: float
    tokamak: Tokamak
    logger: Logger


class NicknameSetting(ABC):
    """
    Abstract base class for nickname settings to resolve tree names for MDSPlus data.

    Methods
    -------
    get_tree_name(params)
        Resolve the tree name based on the given params.
    """

    def get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Resolve the MDSPlus tree name using the provided params.

        Parameters
        ----------
        params : NicknameSettingParams
            The parameters needed to determine the tree name.

        Returns
        -------
        str
            The resolved MDSPlus tree name.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_tree_name(params)

    @abstractmethod
    def _get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Abstract method for resolving the MDSPlus tree name.

        Parameters
        ----------
        params : NicknameSettingParams
            The parameters required to retrieve the tree name.

        Returns
        -------
        str
            The resolved MDSPlus tree name.
        """


class NicknameSettingDict(NicknameSetting):
    """
    A utility class that resolves nicknames using a dictionary of tokamak-nickname
    settings.

    Parameters
    ----------
    nickname_setting_dict : dict[Tokamak, NicknameSettingType]
        A dictionary mapping tokamaks to nickname settings, e.g. `{'cmod': 'efit'}`.
        Any other option passable as a value in the `nickname_setting_dict` dictionary parameter
        in `RetrievalSettings` may be used.
    """

    def __init__(self, nickname_setting_dict: Dict[Tokamak, NicknameSettingType]):
        """
        Initialize with a dictionary of tokamak-nickname mappings.

        Parameters
        ----------
        nickname_setting_dict : dict
            Dictionary of tokamak to nickname mappings.
        """
        resolved_nickname_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_nickname_setting(
                individual_setting
            )
            for tokamak, individual_setting in nickname_setting_dict.items()
        }
        self.resolved_nickname_setting_dict = resolved_nickname_setting_dict

    def _get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Get the tree name based on the resolved nickname setting.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.

        Returns
        -------
        str
            The resolved tree name.
        """
        chosen_setting = self.resolved_nickname_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_tree_name(params)
        raise NotImplementedError(
            f"{self.__class__.__name__} is not implemented for tokamak {params.tokamak}."
        )


class StaticNicknameSetting(NicknameSetting):
    """
    Static nickname setting that always returns a fixed tree name.

    Parameters
    ----------
    tree_name : str
        The tree name to be returned.
    """

    def __init__(self, tree_name: str):
        """
        Initialize StaticNicknameSetting with a fixed tree name.

        Parameters
        ----------
        tree_name : str
            The fixed tree name.
        """
        self.tree_name = tree_name

    def _get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Return the static tree name.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.

        Returns
        -------
        str
            The fixed tree name.
        """
        return self.tree_name


class DefaultNicknameSetting(NicknameSetting):
    """
    Nickname setting to resolve the '_efit_tree' nickname to the default EFIT tree.
    """

    def __init__(self):
        """
        Initialize with default tree names for specific tokamaks.
        """
        self.tokamak_overrides = {
            Tokamak.CMOD: StaticNicknameSetting("analysis")._get_tree_name,
            Tokamak.D3D: StaticNicknameSetting("efit01")._get_tree_name,
        }

    def _get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Raise error if tokamak override is not found.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} is not implemented for tokamak {params.tokamak}."
        )


class DisruptionNicknameSetting(NicknameSetting):
    """
    Nickname setting to resolve the '_efit_tree' nickname to the disruption EFIT tree.
    """

    def __init__(self):
        """
        Initialize with tokamak-specific overrides for disruption trees.
        """
        self.tokamak_overrides = {
            Tokamak.CMOD: self._cmod_nickname,
            Tokamak.D3D: self._d3d_nickname,
        }

    def _d3d_nickname(self, params: NicknameSettingParams) -> str:
        """
        Get the disruption EFIT tree name for D3D.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.

        Returns
        -------
        str
            The resolved EFIT tree name.
        """
        if params.disruption_time is None:
            # TODO: some DIII-D shots have a disruption efit tree, but no disruption time.
            return DefaultNicknameSetting().get_tree_name(params)
        efit_trees = params.database.query(
            "select tree from code_rundb.dbo.plasmas where "
            f"shot = {params.shot_id} and runtag = 'DIS' and deleted = 0 order by idx",
            use_pandas=False,
        )
        if len(efit_trees) == 0:
            return DefaultNicknameSetting().get_tree_name(params)
        efit_tree = efit_trees[-1][0]
        return efit_tree

    def _cmod_nickname(self, params: NicknameSettingParams) -> str:
        """
        Get the disruption EFIT tree name for CMOD.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.

        Returns
        -------
        str
            The resolved EFIT tree name.
        """
        if params.disruption_time is None:
            return DefaultNicknameSetting().get_tree_name(params)
        return "efit18"

    def _get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Raise error if tokamak override is not found.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.
        """
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
    Resolve a nickname setting to its corresponding NicknameSetting object.

    Parameters
    ----------
    nickname_setting : NicknameSettingType
        The nickname setting, which can be a string, a dictionary, or a NicknameSetting instance.

    Returns
    -------
    NicknameSetting
        The resolved NicknameSetting object.
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
