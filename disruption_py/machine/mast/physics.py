#!/usr/bin/env python3

"""
Physics methods for MAST.
"""

import xarray as xr
from s3fs import S3Map
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.inout.s3 import S3Connection
from disruption_py.machine.tokamak import Tokamak


class MastPhysicsMethods:
    """
    MAST class to serve as a method holder.
    """

    @staticmethod
    @physics_method(
        columns=["ip"],
        tokamak=Tokamak.MAST,
    )
    def get_ip_parameters(params: PhysicsMethodParams):
        """Get Ip parameters"""
        conn: S3Connection = params.mds_conn.conn
        file_name = f"s3://mast/level2/shots/{params.shot_id}.zarr"
        mapper = conn.get_mapper(file_name)
        ds = xr.open_zarr(mapper, group="summary")
        ip = ds["ip"].values
        return {"ip": ip}
