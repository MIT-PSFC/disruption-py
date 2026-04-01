#!/usr/bin/env python3

"""
Abstract base classes for data connections.

DataConnection: per-shot data access (get_data, get_data_with_dims, get_dims).
ProcessConnection: per-process factory that creates DataConnection instances.
"""

from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np

from disruption_py.machine.tokamak import Tokamak


class DataConnection(ABC):
    """Per-shot data access interface.

    Each instance is bound to a single shot. Implementations must provide
    get_data, get_data_with_dims, get_dims, and cleanup. The reconnect
    method is optional (default no-op).
    """

    @property
    @abstractmethod
    def shot_id(self) -> int:
        """The shot ID this connection is bound to."""

    @abstractmethod
    def get_data(self, path: str, group: str = None, **kwargs) -> np.ndarray:
        """Get data at path.

        Parameters
        ----------
        path : str
            Data path (node path for MDSplus, variable name for Xarray).
        group : str, optional
            Container name (tree for MDSplus, group for Xarray).
        **kwargs
            Backend-specific options.

        Returns
        -------
        np.ndarray
        """

    @abstractmethod
    def get_data_with_dims(
        self,
        path: str,
        group: str = None,
        dim_nums: List = None,
        **kwargs,
    ) -> Tuple:
        """Get data and dimension arrays.

        Parameters
        ----------
        path : str
            Data path.
        group : str, optional
            Container name.
        dim_nums : List, optional
            Dimension indices to retrieve. Default [0].
        **kwargs
            Backend-specific options.

        Returns
        -------
        Tuple
            (data, dim0, dim1, ...) as numpy arrays.
        """

    @abstractmethod
    def get_dims(
        self,
        path: str,
        group: str = None,
        dim_nums: List = None,
        **kwargs,
    ) -> Tuple:
        """Get only dimension arrays.

        Parameters
        ----------
        path : str
            Data path.
        group : str, optional
            Container name.
        dim_nums : List, optional
            Dimension indices to retrieve. Default [0].
        **kwargs
            Backend-specific options.

        Returns
        -------
        Tuple
            Requested dimensions.
        """

    @abstractmethod
    def cleanup(self) -> None:
        """Release resources for this shot."""

    def reconnect(self) -> None:
        """Reconnect after error. Default no-op.

        Xarray opens a new DataTree per shot so there is nothing to
        reconnect. MDSplus overrides this to call conn.reconnect().
        """


class ProcessConnection(ABC):
    """Per-process factory that creates DataConnection instances."""

    @abstractmethod
    def get_shot_connection(self, shot_id: int) -> DataConnection:
        """Create a per-shot DataConnection for the given shot."""

    @classmethod
    @abstractmethod
    def from_config(cls, tokamak: Tokamak) -> "ProcessConnection":
        """Create a ProcessConnection from tokamak configuration."""
