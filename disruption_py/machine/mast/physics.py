#!/usr/bin/env python3

"""
Physics methods for MAST.
"""

import xarray as xr
from disruption_py.core.utils.math import interp1
from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.inout.xarray import XarrayConnection
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
        conn: XarrayConnection = params.mds_conn
        file_name = conn.get_shot_file_path(params.shot_id)
        ds = xr.open_zarr(file_name, group="summary")
        ip = ds["ip"].values
        magtime = ds["time"].values

        times = params.times
        ip = interp1(magtime, ip, times)
        return {"ip": ip}
