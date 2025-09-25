#!/usr/bin/env python3

"""
Module for retrieving and calculating data for HBTEP physics methods.
"""

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak


class HbtepPhysicsMethods:
    """
    This class provides methods to retrieve and calculate physics-related data
    for HBTEP.
    """

    @staticmethod
    @physics_method(columns=["ip"], tokamak=Tokamak.HBTEP)
    def get_ip(params: PhysicsMethodParams):
        """
        Get the plasma current
        """
        ip, t_ip = params.mds_conn.get_data_with_dims(
            r"\TOP.SENSORS.ROGOWSKIS:IP", tree_name="hbtep2"
        )
        ip = interp1(t_ip, ip, params.times, "linear")
        return {"ip": ip}
