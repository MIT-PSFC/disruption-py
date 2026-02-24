#!/usr/bin/env python3

"""
Module for connecting to Xarray data stores.
"""

import threading
from typing import List, Tuple

import numpy as np
import xarray as xr
from loguru import logger

from disruption_py.config import config
from disruption_py.core.utils.misc import shot_msg
from disruption_py.inout.base import DataConnection, ProcessConnection
from disruption_py.machine.tokamak import Tokamak


class ProcessXarrayConnection(ProcessConnection):
    """
    Process-level Xarray connection.

    Holds configuration for the Xarray store and creates per-shot
    XarrayDataConnection instances.
    """

    def __init__(
        self,
        file_path: str,
        file_ext: str = "zarr",
        endpoint_url: str | None = None,
    ):
        if file_path is None:
            return

        self.endpoint_url = endpoint_url
        self.file_path = file_path
        self.file_ext = file_ext

        logger.debug(
            "PID #{pid} | Connecting to Xarray store: {server}",
            server=endpoint_url,
            pid=threading.get_native_id(),
        )

    @property
    def folder_path(self):
        """Get full folder path including endpoint URL if provided."""
        if self.endpoint_url is not None:
            return f"{self.endpoint_url}/{self.file_path}"

        return self.file_path

    @classmethod
    def from_config(cls, tokamak: Tokamak) -> "ProcessXarrayConnection":
        """
        Create instance of the connection based on the file path, file extension,
        and endpoint URL from the configuration.
        """
        params = config(tokamak).inout.xarray
        return ProcessXarrayConnection(**params)

    def get_shot_file_path(self, shot_id: int):
        """Get file path for individual shot."""
        return f"{self.folder_path}/{shot_id}.{self.file_ext}"

    def get_shot_connection(self, shot_id: int) -> "XarrayDataConnection":
        """Get per-shot data connection for individual shot."""
        file_path = self.get_shot_file_path(shot_id)
        engine = "zarr" if self.file_ext == "zarr" else "netcdf4"

        logger.trace(
            shot_msg("Opening data tree: {file_path}"),
            shot=shot_id,
            file_path=file_path,
        )

        data_tree = xr.open_datatree(
            file_path, engine=engine, chunks=None, create_default_indexes=False
        )
        return XarrayDataConnection(shot_id, data_tree)


class XarrayDataConnection(DataConnection):
    """
    Per-shot Xarray data connection wrapping a DataTree.
    """

    def __init__(self, shot_id: int, data_tree: xr.DataTree):
        self._shot_id = shot_id
        self.data_tree = data_tree

    @property
    def shot_id(self) -> int:
        """The shot ID this connection is bound to."""
        return self._shot_id

    def _resolve_path(self, path: str, group: str = None) -> str:
        """Prepend group to path if group is provided and path has no '/'."""
        if group is not None and "/" not in path:
            return f"{group}/{path}"
        return path

    def get_data(
        self, path: str, group: str = None, return_xarray: bool = False, **kwargs
    ) -> np.ndarray:
        """Get data from the connection.

        Parameters
        ----------
        path : str
            Variable path, e.g. "summary/ip" or just "ip" if group is provided.
        group : str, optional
            Group prefix. If provided and path has no "/", prepends group.
        return_xarray : bool, optional
            If True, return the raw xarray DataArray instead of numpy values.
        **kwargs
            Backend-specific options (ignored).

        Returns
        -------
        np.ndarray or xr.DataArray
        """
        path = self._resolve_path(path, group)
        logger.trace(shot_msg("Getting data: {path}"), shot=self._shot_id, path=path)

        try:
            item = self.data_tree[path]

            if return_xarray:
                return item

            return item.values
        except KeyError:
            logger.warning(
                shot_msg("Variable not found: {path}"),
                path=path,
                shot=self._shot_id,
            )

        if return_xarray:
            return None

        return np.array([np.nan])

    def get_data_with_dims(
        self,
        path: str,
        group: str = None,
        dim_nums: List = None,
        **kwargs,
    ) -> Tuple:
        """Get data and dimension arrays from the DataTree.

        Parameters
        ----------
        path : str
            Variable path.
        group : str, optional
            Group prefix.
        dim_nums : List, optional
            Dimension indices to retrieve. Default [0].
        **kwargs
            Backend-specific options (ignored).

        Returns
        -------
        Tuple
            (data, dim0, dim1, ...) as numpy arrays.
        """
        path = self._resolve_path(path, group)
        dim_nums = dim_nums or [0]
        logger.trace(
            shot_msg("Getting data and dims: {path}"), shot=self._shot_id, path=path
        )

        item = self.data_tree[path]
        data = item.values
        dim_names = list(item.dims)
        dims = []
        for d in dim_nums:
            dim_name = dim_names[d]
            dims.append(item.coords[dim_name].values)

        return data, *dims

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
            Variable path.
        group : str, optional
            Group prefix.
        dim_nums : List, optional
            Dimension indices to retrieve. Default [0].
        **kwargs
            Backend-specific options (ignored).

        Returns
        -------
        Tuple
            Requested dimensions.
        """
        result = self.get_data_with_dims(path, group=group, dim_nums=dim_nums, **kwargs)
        return result[1:]

    def cleanup(self) -> None:
        """Close the DataTree."""
        if self.data_tree is not None:
            self.data_tree.close()
            self.data_tree = None


# Backward compatibility alias
XarrayConnection = ProcessXarrayConnection
