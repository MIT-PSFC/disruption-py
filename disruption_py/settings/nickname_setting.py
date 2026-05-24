#!/usr/bin/env python3

"""
This module provides classes for handling nickname settings, which are used to
resolve MDSplus tree names for various tokamaks.
"""

import os
import sys
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Union

import pandas as pd
from loguru import logger

from disruption_py.config import config
from disruption_py.core.physics_method.errors import FetchDataError
from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.core.utils.misc import get_temporary_folder
from disruption_py.inout.base import DataConnection
from disruption_py.inout.mds import mdsExceptions
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak

NicknameSettingType = Union[
    "NicknameSetting",
    str,
    List[Union[str, "NicknameSetting"]],
    Dict[Tokamak, "NicknameSettingType"],
]


@dataclass
class NicknameSettingParams:
    """
    Parameters passed to nickname trees for resolving tree nicknames.

    Attributes
    ----------
    shot_id : int
        The shot ID for which to resolve nicknames.
    data_conn : DataConnection
        Data connection for the shot.
    database : ShotDatabase
        Database connection for querying tokamak shot data.
    disruption_time : float
        The time of the disruption in seconds.
    tokamak : Tokamak
        The tokamak for which results are being processed.
    """

    shot_id: int
    data_conn: DataConnection
    database: ShotDatabase
    disruption_time: float
    tokamak: Tokamak


def _efit_df_to_dict(df: pd.DataFrame) -> dict[int, str]:
    """
    Build shot/tree dict from an indexed DataFrame.
    When multiple rows exist for the same shot, select the largest idx.
    """
    df = (
        df.reset_index()
        .sort_values("idx", ascending=False)
        .drop_duplicates("shot", keep="first")
        .set_index("shot")["tree"]
    )
    return df.to_dict()


class NicknameSetting(ABC):
    """
    Abstract base class for nickname settings to resolve tree names for MDSPlus data.

    Methods
    -------
    get_tree_name(params)
        Resolve the tree name based on the given params.
    prefetch_db(database, tokamak)
        Pre-fetch EFIT trees for all shots, if available.
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

    def prefetch_db(self, database: ShotDatabase, tokamak: Tokamak) -> None:
        """
        Pre-fetch EFIT trees for all shots, if available.

        Parameters
        ----------
        database : ShotDatabase
            Database connection for querying tokamak shot data.
        tokamak : Tokamak
            The tokamak for which results are being processed.
        """

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

    def prefetch_db(self, database: ShotDatabase, tokamak: Tokamak) -> None:
        """
        Map pre-fetching to the tokamak-specific nickname setting.

        Parameters
        ----------
        database : ShotDatabase
            Database connection for querying tokamak shot data.
        tokamak : Tokamak
            The tokamak for which results are being processed.
        """
        setting = self.resolved_nickname_setting_dict.get(tokamak)
        if setting is not None:
            setting.prefetch_db(database, tokamak)

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


class NicknameSettingList(NicknameSetting):
    """
    Cascade nickname setting. Tries each item in order via ``open_tree``,
    returning the first that opens successfully. Falls back to the next on
    ``mdsExceptions.TreeFOPENR``.

    Parameters
    ----------
    items : list
        Cascade of tree-name strings or NicknameSetting instances. Tried in
        order; the first item that opens wins.
    """

    def __init__(self, items: list):
        """
        Initialize with a cascade of items to try in order.

        Parameters
        ----------
        items : list
            Non-empty cascade. Each item is a tree-name string or a
            NicknameSetting instance.
        """
        if not items:
            raise ValueError("NicknameSettingList requires at least one item.")
        self.items = items
        self.resolved_items = [self._resolve_item(item) for item in items]

    @staticmethod
    def _resolve_item(item):
        """Validate a cascade item: must be a string or NicknameSetting."""
        if isinstance(item, NicknameSetting):
            return item
        if isinstance(item, str):
            return item
        raise ValueError(
            "NicknameSettingList items must be str or NicknameSetting, "
            f"got {type(item).__name__}: {item!r}"
        )

    def _get_tree_name(self, params: NicknameSettingParams) -> str:
        """
        Try each cascade item via ``open_tree``; return the first that opens.

        Parameters
        ----------
        params : NicknameSettingParams
            Parameters needed to determine the nickname.

        Returns
        -------
        str
            The tree name of the first cascade item that opened successfully.

        Raises
        ------
        FetchDataError
            If no item in the cascade could be opened.
        """
        attempts: List[str] = []
        for item in self.resolved_items:
            candidate = item if isinstance(item, str) else type(item).__name__
            try:
                if isinstance(item, NicknameSetting):
                    candidate = item.get_tree_name(params)
                params.data_conn.open_tree(candidate)
            except (mdsExceptions.TreeFOPENR, FetchDataError):
                # TreeFOPENR = tree missing; FetchDataError = nested cascade exhausted.
                attempts.append(candidate)
                continue
            if attempts:
                logger.verbose(
                    "Nickname cascade for shot {shot}: selected '{name}' "
                    "after failed: {prior}",
                    shot=params.shot_id,
                    name=candidate,
                    prior=", ".join(attempts),
                )
            return candidate
        raise FetchDataError(
            f"No tree in nickname cascade could be opened for shot "
            f"{params.shot_id}. Tried: {attempts}"
        )


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
            Tokamak.EAST: StaticNicknameSetting("efit_east")._get_tree_name,
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

    EAST does not have a disruption EFIT tree, so we'll use the default "efit_east" tree.
    """

    def __init__(self):
        """
        Initialize with tokamak-specific overrides for disruption trees.
        """
        self.tokamak_overrides = {
            Tokamak.CMOD: self._cmod_nickname,
            Tokamak.D3D: self._d3d_nickname,
            Tokamak.EAST: StaticNicknameSetting("efit_east")._get_tree_name,
        }
        self._d3d_efit_trees: dict | None = None

    def prefetch_db(self, database: ShotDatabase, tokamak: Tokamak) -> None:
        """
        Pre-fetch EFIT trees for all DIII-D shots in a single query.
        Results are cached to a CSV file in the daily temporary folder.

        Parameters
        ----------
        database : ShotDatabase
            Database connection for querying tokamak shot data.
        tokamak : Tokamak
            The tokamak for which results are being processed.
        """
        if tokamak != Tokamak.D3D:
            return

        runtag = config(tokamak).efit.runtag
        if "pytest" in sys.modules:
            runtag = "DIS"

        # write to unique runid folder, cache to parent
        tmp = os.path.join(get_temporary_folder(), f"rundb_{runtag}.csv")
        csv = os.path.join(
            os.path.dirname(get_temporary_folder()), f"rundb_{runtag}.csv"
        )

        # cache hit
        if os.path.isfile(csv):
            logger.debug(
                "Loading EFIT '{runtag}' trees from cache: {csv}",
                runtag=runtag,
                csv=csv,
            )
            df = pd.read_csv(csv, index_col="idx")
            self._d3d_efit_trees = _efit_df_to_dict(df)
            return

        # cache miss
        logger.debug("Fetching EFIT '{runtag}' trees from DB...", runtag=runtag)
        took = -time.time()
        df = database.query(
            "select idx, shot, tree "
            "from code_rundb.dbo.plasmas "
            f"where runtag = '{runtag}' and deleted = 0 "
            "order by idx"
        )
        took += time.time()

        # empty
        if df is None or df.empty:
            logger.warning("Fetched empty EFIT '{runtag}' trees.", runtag=runtag)
            self._d3d_efit_trees = {}
            return

        # log results
        logger.debug(
            "Fetched EFIT '{runtag}' trees in {sec:.3f}s: {tot:,} rows, {unique:,} unique shots.",
            runtag=runtag,
            sec=took,
            tot=len(df),
            unique=df["shot"].nunique(),
        )

        # store
        df = df.set_index("idx")
        self._d3d_efit_trees = _efit_df_to_dict(df)

        # write cache
        logger.debug("Caching EFIT '{runtag}' trees: {csv}", runtag=runtag, csv=csv)
        df.to_csv(tmp)
        if os.path.exists(csv):
            logger.warning("Removing stale EFIT cache: {csv}", csv=csv)
            os.remove(csv)
        os.link(tmp, csv)

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

        # all shots have "efit21" / DISPY tree, but
        # only disruptive shots have "efit18" / DIS tree, thus
        # non-disruptive shots should use the default tree

        runtag = config(params.tokamak).efit.runtag
        if "pytest" in sys.modules:
            runtag = "DIS"
        if runtag == "DIS" and params.disruption_time is None:
            return DefaultNicknameSetting().get_tree_name(params)

        # prefetch if possible
        if self._d3d_efit_trees is not None:
            tree = self._d3d_efit_trees.get(params.shot_id)
            if tree is not None:
                return tree
            return DefaultNicknameSetting().get_tree_name(params)

        # fallback single query
        efit_trees = params.database.query(
            "select tree from code_rundb.dbo.plasmas where "
            f"shot = {params.shot_id} and runtag = '{runtag}' and deleted = 0 order by idx",
            use_pandas=False,
        )
        if len(efit_trees) == 0:
            return DefaultNicknameSetting().get_tree_name(params)
        return efit_trees[-1][0]

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

        # all shots have efit21 tree, but
        # only disruptive shots have efit18 tree, thus
        # non-disruptive shots should use the default tree

        tree = config(params.tokamak).efit.tree
        if "pytest" in sys.modules:
            tree = "efit18"
        return NicknameSettingList([tree, "analysis"]).get_tree_name(params)

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
        The nickname setting: a NicknameSetting instance, a string (single tree
        name, registered key like "disruption", or comma-separated cascade), a
        list (cascade), or a dictionary (per-tokamak dispatch).

    Returns
    -------
    NicknameSetting
        The resolved NicknameSetting object.
    """
    if isinstance(nickname_setting, NicknameSetting):
        return nickname_setting
    if isinstance(nickname_setting, list):
        return NicknameSettingList(nickname_setting)
    if isinstance(nickname_setting, dict):
        return NicknameSettingDict(nickname_setting)
    if isinstance(nickname_setting, str):
        if nickname_setting in _nickname_setting_mappings:
            return _nickname_setting_mappings[nickname_setting]
        if "," in nickname_setting:
            tokens = [s.strip() for s in nickname_setting.split(",") if s.strip()]
            # Registered keys map to their class; other tokens are literal tree names.
            items = [
                _nickname_setting_mappings[s] if s in _nickname_setting_mappings else s
                for s in tokens
            ]
            return NicknameSettingList(items)
        return StaticNicknameSetting(nickname_setting)

    raise ValueError(f"Invalid nickname setting type {type(nickname_setting)}.")
