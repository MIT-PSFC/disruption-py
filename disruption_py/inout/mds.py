#!/usr/bin/env python3

"""
Module for managing connections to MDSplus.
"""

from typing import Any, Callable, Dict, List, Tuple

import MDSplus
import numpy as np
from loguru import logger

from disruption_py.config import config
from disruption_py.core.utils.misc import make_hashable, safe_cast
from disruption_py.core.utils.shared_instance import SharedInstance
from disruption_py.machine.tokamak import Tokamak


class ProcessMDSConnection:
    """
    Abstract class for connecting to MDSplus.

    Ensure that a single MDSPlus connection is used by each process for all shots
    retrieved by that process.
    """

    def __init__(self, conn_string: str):
        self.conn = None
        if conn_string is None:
            return
        # pylint: disable=no-member
        self.conn = MDSplus.Connection(conn_string)
        try:
            self.conn.get("shorten_path()")
        except MDSplus.mdsExceptions.TdiUNKNOWN_VAR:
            logger.debug("MDSplus does not support the `shorten_path()` method.")

    @classmethod
    def from_config(cls, tokamak: Tokamak):
        """
        Create instance of the MDS connection based on the connection string
        from the configuration.
        """
        return SharedInstance(ProcessMDSConnection).get_instance(
            config(tokamak).inout.mds.mdsplus_connection_string
        )

    def get_shot_connection(self, shot_id: int):
        """Get MDSPlus Connection wrapper for individual shot."""
        return MDSConnection(self.conn, shot_id)


def _better_mds_exceptions(func):
    """
    Decorator to catch MDSplus exceptions and recreate them with a more descriptive
    error message which includes the tree and the node path.
    """

    def wrapper(*args, **kwargs):
        self = args[0]
        if len(args) > 1:
            path = args[1]
        else:
            path = kwargs.get("path", None)
        try:
            return func(*args, **kwargs)
        except MDSplus.mdsExceptions.TreeFOPENR:
            nick = kwargs.get("tree_name", None)
            tree = self.tree_name(nick)
            if nick == tree:
                nick = ""
            raise MDSplus.mdsExceptions.TreeFOPENR(
                "Tree not found. "
                + (f"Nick: {nick}, " if nick else "")
                + f"Tree: {tree}"
            ) from None
        except MDSplus.mdsExceptions.TreeNNF:
            raise MDSplus.mdsExceptions.TreeNNF(
                f"Node not found. Tree: {self.last_open_tree}, Node: {path}"
            ) from None
        except MDSplus.mdsExceptions.TreeNODATA:
            raise MDSplus.mdsExceptions.TreeNODATA(
                f"No data available. Tree: {self.last_open_tree}, Node: {path}"
            ) from None

    return wrapper


class MDSConnection:
    """
    Wrapper class for MDSPlus Connection class used for handling individual shots.
    """

    def __init__(
        self, conn: MDSplus.Connection, shot_id: int  # pylint: disable=no-member
    ):
        self.conn = conn
        self.shot_id = shot_id
        self.tree_nickname_funcs = {}
        self.tree_nicknames = {}
        self.last_open_tree = None
        self.data_cache = {}

    @_better_mds_exceptions
    def open_tree(self, tree_name: str):
        """
        Open the specified _name.

        If the specified tree_name is a nickname for a tree_name, will open the tree
        that it is a nickname for.
        """
        logger.trace("Opening tree: {tree_name}", tree_name=tree_name)
        if (
            tree_name not in self.tree_nicknames
            and tree_name in self.tree_nickname_funcs
        ):
            self.tree_nicknames[tree_name] = self.tree_nickname_funcs[tree_name]()

        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]

        if self.last_open_tree != tree_name:
            self.conn.openTree(tree_name, self.shot_id)

        self.last_open_tree = tree_name

    def cleanup(self):
        """
        Close all open trees
        """
        self.last_open_tree = None

    @_better_mds_exceptions
    def get(self, expression: str, arguments: Any = None, tree_name: str = None) -> Any:
        """
        Evaluate the specified expression.

        The expression is passed as string argument, but may contain optional arguments.
        These arguments are then passed as an array of Data objects.

        Parameters
        ----------
        expression : str
            MDSplus TDI expression. Please see MDSplus documentation for more
            information.
        arguments : Any, optional
            Arguments for MDSplus TDI Expression. Please see MDSplus documentation
            for more information. Default None.
        tree_name : str, optional

        Returns
        -------
        Any
            Result of evaluating TDI expression from MDSplus.
        """
        # Try to retrieve evaluated expression from cache
        cache_key = make_hashable([expression, tree_name, arguments])
        eval_expr = None

        if cache_key in self.data_cache:
            eval_expr = self.data_cache[cache_key]

        # Retrieve evaluated expression from MDSplus
        if eval_expr is None:
            if tree_name is not None:
                self.open_tree(tree_name)
            if arguments is None:
                eval_expr = self.conn.get(expression)
            else:
                eval_expr = self.conn.get(expression, arguments)

            # Cache evaluated expression
            if cache_key not in self.data_cache:
                self.data_cache[cache_key] = {}
            self.data_cache[cache_key] = eval_expr

        return eval_expr

    # Convenience methods

    @_better_mds_exceptions
    def get_data(
        self,
        path: str,
        tree_name: str = None,
        astype: str = None,
        arguments: Any = None,
    ) -> np.ndarray:
        """
        Get data for record at specified path.

        Parameters
        ----------
        path : str
            MDSplus path to record.
        tree_name : str, optional
            The name of the tree that must be open for retrieval.
        astype : str, optional
            The data type for explicit casting.
        arguments : Any, optional
            Arguments for MDSplus TDI Expression. Default None.
            Please see MDSplus documentation for more information.

        Returns
        -------
        np.ndarray
            Returns the node data.
        """
        # Try to retrieve data from cache
        cache_key = make_hashable([path, tree_name, arguments])

        data = None
        if cache_key in self.data_cache:
            data = self.data_cache[cache_key].get("data")

        if data is None:
            # Retrieve data from MDSplus
            if tree_name is not None:
                self.open_tree(tree_name)

            logger.trace("Getting data: {path}", path=path)
            data = self.conn.get(path, arguments).data()

            # Cache data
            if cache_key not in self.data_cache:
                self.data_cache[cache_key] = {}
            self.data_cache[cache_key]["data"] = data

        if astype:
            data = safe_cast(data, astype)

        return data

    @_better_mds_exceptions
    def get_data_with_dims(
        self,
        path: str,
        tree_name: str = None,
        dim_nums: List = None,
        astype: str = None,
        cast_all: bool = False,
    ) -> Tuple:
        """
        Get data and dimension(s) for record at specified path.

        Parameters
        ----------
        path : str
            MDSplus path to record.
        tree_name : str, optional
            The name of the tree that must be open for retrieval.
        dim_nums : List, optional
            A list of dimensions that should have their size retrieved. Default [0].
        astype : str, optional
            The data type for explicit casting.
        cast_all : bool, optional. Default False.
            Whether to cast both data and dims, or only data.

        Returns
        -------
        Tuple
            Returns the node data, followed by the requested dimensions.
        """

        dim_nums = dim_nums or [0]

        if tree_name is not None:
            self.open_tree(tree_name)

        logger.trace("Getting data and dims: {path}", path=path)
        data = self.get_data(path=path, tree_name=tree_name, astype=astype)
        dims = self.get_dims(
            path=path,
            tree_name=tree_name,
            dim_nums=dim_nums,
            astype=astype if cast_all else None,
        )

        return data, *dims

    @_better_mds_exceptions
    def get_dims(
        self,
        path: str,
        tree_name: str = None,
        dim_nums: List = None,
        astype: str = None,
    ) -> Tuple:
        """
        Get the specified dimensions (e.g. timebase) for record at specified path.

        Parameters
        ----------
        path : str
            MDSplus path to record.
        tree_name : str, optional
            The name of the tree that must be open for retrieval.
        dim_nums : List, optional
            A list of dimensions that should have their size retrieved. Default [0].
        astype : str, optional
            The data type for explicit casting.

        Returns
        -------
        Tuple
            Returns the requested dimensions as a tuple.
        """
        dim_nums = dim_nums or [0]
        dims = [None for _ in dim_nums]

        # Try to retrieve dimensions of data from cache.
        cache_key = make_hashable([path, tree_name])
        if cache_key in self.data_cache:
            cached: dict = self.data_cache[cache_key]
            for i, dim_num in enumerate(dim_nums):
                if dim_num in cached:
                    dims[i] = cached[dim_num]

        # Retrieve from MDSplus
        retrieved_from_cache = any(d is not None for d in dims)
        if not retrieved_from_cache:
            if tree_name is not None:
                self.open_tree(tree_name)

            logger.trace("Getting dims: {path}", path=path)
            dims = [self.conn.get(f"dim_of({path},{d})").data() for d in dim_nums]

            # Cache dims
            if cache_key not in self.data_cache:
                self.data_cache[cache_key] = {}
            for i, dim_num in enumerate(dim_nums):
                self.data_cache[cache_key][dim_num] = dims[i]
        if astype:
            dims = [safe_cast(dim, astype) for dim in dims]

        return dims

    # nicknames

    def add_tree_nickname_funcs(self, tree_nickname_funcs: Dict[str, Callable]):
        """
        Add tree nickname functions to the connection.

        Required because some tree nickname functions require the connection to exist.
        """
        self.tree_nickname_funcs.update(tree_nickname_funcs)

    def get_tree_name_of_nickname(self, nickname: str):
        """
        Get the tree name that the nickname has been set to or None if the nickname
        was not set.
        """
        if nickname not in self.tree_nicknames and nickname in self.tree_nickname_funcs:
            self.tree_nicknames[nickname] = self.tree_nickname_funcs[nickname]()

        return self.tree_nicknames.get(nickname, None)

    def tree_name(self, for_name: str) -> str:
        """
        The tree name for for_name, whether it is a nickname or tree name itself
        """
        return self.get_tree_name_of_nickname(for_name) or for_name
