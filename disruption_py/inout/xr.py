#!/usr/bin/env python3

"""
Module for connecting to Xarray data stores.
"""

import threading

import numpy as np
import xarray as xr
from loguru import logger

from disruption_py.config import config
from disruption_py.core.utils.misc import shot_msg
from disruption_py.machine.tokamak import Tokamak


class XarrayConnection:
    """
    Class for connecting to Xarray store.
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
        self.data_tree: xr.DataTree | None = None

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
    def from_config(cls, tokamak: Tokamak):
        """
        Create instance of the connection based on the file path, file extension, and endpoint URL
        from the configuration.
        """
        params = config(tokamak).inout.xarray
        return XarrayConnection(**params)

    def get_shot_connection(self, shot_id: int):
        """Get connection to xarray store for individual shot."""
        file_path = self.get_shot_file_path(shot_id)
        engine = "zarr" if self.file_ext == "zarr" else "netcdf4"

        logger.trace(
            shot_msg("Opening data tree: {file_path}"),
            shot=shot_id,
            file_path=file_path,
        )

        self.data_tree = xr.open_datatree(
            file_path, engine=engine, chunks=None, create_default_indexes=False
        )
        return self

    def get_shot_file_path(self, shot_id: int):
        """Get file path for individual shot."""
        file_path = f"{self.folder_path}/{shot_id}.{self.file_ext}"
        return file_path

    def get_data(self, shot_id: int, path: str, return_xarray: bool = False):
        """Get data from the connection."""
        if self.data_tree is None:
            self.get_shot_connection(shot_id)

        try:
            item = self.data_tree[path]

            if return_xarray:
                return item

            return item.values
        except KeyError:
            logger.warning(
                shot_msg("Variable not found: {path}"),
                path=path,
                shot=shot_id,
            )

        if return_xarray:
            return None

        return np.array([np.nan])

    def cleanup(self):
        """Cleanup the connection."""

    def reconnect(self):
        """Reconnect the connection."""
