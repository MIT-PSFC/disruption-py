#!/usr/bin/env python3

"""
Module for managing connections to MDSplus.
"""

import logging
from typing import Any, Callable, Dict, List, Tuple

import MDSplus
import numpy as np

from disruption_py.config import config
from disruption_py.core.utils.misc import safe_cast
from disruption_py.core.utils.shared_instance import SharedInstance
from disruption_py.machine.tokamak import Tokamak


class ProcessMDSConnection:
    """
    Abstract class for connecting to MDSplus.

    Ensure that a single MDSPlus connection is used by each process for all shots
    retrieved by that process.
    """

    logger = logging.getLogger("disruption_py")

    def __init__(self, conn_string: str):
        self.conn = None
        if conn_string is None:
            return
        # pylint: disable=no-member
        self.conn = MDSplus.Connection(conn_string)
        try:
            self.conn.get("shorten_path()")
        except MDSplus.mdsExceptions.TdiUNKNOWN_VAR:
            self.logger.debug("MDSplus does not support the `shorten_path()` method.")

    @classmethod
    def from_config(cls, tokamak: Tokamak):
        """
        Create instance of the MDS connection based on the connection string
        from the configuration.
        """
        return SharedInstance(ProcessMDSConnection).get_instance(
            config(tokamak).MDSPLUS_CONNECTION_STRING
        )

    def get_shot_connection(self, shot_id: int):
        """Get MDSPlus Connection wrapper for individual shot."""
        return MDSConnection(self.conn, shot_id)


class MDSConnection:
    """
    Wrapper class for MDSPlus Connection class used for handling individual shots.
    """

    logger = logging.getLogger("disruption_py")

    def __init__(
        self, conn: MDSplus.Connection, shot_id: int  # pylint: disable=no-member
    ):
        self.conn = conn
        self.shot_id = shot_id
        self.tree_nickname_funcs = {}
        self.tree_nicknames = {}
        self.last_open_tree = None

    def open_tree(self, tree_name: str):
        """
        Open the specified _name.

        If the specified tree_name is a nickname for a tree_name, will open the tree
        that it is a nickname for.
        """
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
        if tree_name is not None:
            self.open_tree(tree_name)
        if arguments is None:
            return self.conn.get(expression)
        return self.conn.get(expression, arguments)

    # Convenience methods

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

        if tree_name is not None:
            self.open_tree(tree_name)

        data = self.conn.get("_sig=" + path, arguments).data()
        if astype:
            data = safe_cast(data, astype)

        return data

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

        data = self.conn.get("_sig=" + path).data()
        dims = [self.conn.get(f"dim_of(_sig,{dim_num})").data() for dim_num in dim_nums]

        if astype:
            data = safe_cast(data, astype)
            if cast_all:
                dims = [safe_cast(dim, astype) for dim in dims]

        return data, *dims

    def get_dims(
        self,
        path: str,
        tree_name: str = None,
        dim_nums: List = None,
        astype: str = None,
    ) -> Tuple:
        """
        Get the specified dimensions for record at specified path.

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

        if tree_name is not None:
            self.open_tree(tree_name)

        dims = [self.conn.get(f"dim_of({path},{d})").data() for d in dim_nums]

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
