import threading
import numpy as np
import xarray as xr
from typing import Optional
from loguru import logger
from zarr import group

from disruption_py.machine.tokamak import Tokamak
from disruption_py.config import config


class XarrayConnection:
    """
    Class for connecting to Xarray store.
    """

    def __init__(
        self,
        file_path: str,
        file_ext: str = "zarr",
        endpoint_url: Optional[str] = None,
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
        # pylint: disable=no-member

    @property
    def folder_path(self):
        if self.endpoint_url is not None:
            return f"{self.endpoint_url}/{self.file_path}"
        else:
            return self.file_path

    @classmethod
    def from_config(cls, tokamak: Tokamak):
        """
        Create instance of the connection based on the file path, file extension, and endpoint URL
        from the configuration.
        """
        params = dict(config(tokamak).inout.xarray)
        return XarrayConnection(**params)

    def get_shot_connection(self, shot_id: int):
        """Get connection."""
        file_path = self.get_shot_file_path(shot_id)
        engine = "zarr" if self.file_ext == "zarr" else "netcdf4"
        self.data_tree = xr.open_datatree(file_path, engine=engine)
        return self

    def get_shot_file_path(self, shot_id: int):
        """Get file path for individual shot."""
        file_path = f"{self.folder_path}/{shot_id}.{self.file_ext}"
        return file_path

    def get_data(self, shot_id: int, path: str):
        """Get data from the connection."""
        if self.data_tree is None:
            self.get_shot_connection(shot_id)

        try:
            return self.data_tree[path].values
        except KeyError:
            logger.warning(
                "PID #{pid} | Variable '{path}' not found for shot {shot_id}.",
                pid=threading.get_native_id(),
                path=path,
                shot_id=shot_id,
            )
            return np.array([np.nan])


    def cleanup(self):
        pass
