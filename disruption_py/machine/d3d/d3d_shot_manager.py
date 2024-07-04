#!/usr/bin/env python3

import traceback

import numpy as np

from disruption_py.io.sql import ShotDatabase
from disruption_py.io.mds import (
    MDSConnection,
    ProcessMDSConnection,
)
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.machine.shot_manager import ShotManager
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak
from disruption_py.core.utils.math import interp1
from disruption_py.core.utils.misc import without_duplicates


class D3DShotManager(ShotManager):

    @property
    def tokamak(self) -> Tokamak:
        return Tokamak.D3D

    @classmethod
    def _modify_times_flattop_timebase(
        cls, physics_method_params: PhysicsMethodParams, **kwargs
    ):
        try:
            (
                ip_prog,
                t_ip_prog,
            ) = physics_method_params.mds_conn.get_data_with_dims(
                f"ptdata('iptipp', {physics_method_params.shot_id})", tree_name="d3d"
            )
            t_ip_prog = t_ip_prog / 1.0e3  # [ms] -> [s]
            polarity = np.unique(
                physics_method_params.mds_conn.get_data(
                    f"ptdata('iptdirect', {physics_method_params.shot_id})",
                    tree_name="d3d",
                )
            )
            if len(polarity) > 1:
                cls.logger.info(
                    f"[Shot {physics_method_params.shot_id}]:Polarity of Ip target is not constant. Using value at first timestep."
                )
                cls.logger.debug(
                    f"[Shot {physics_method_params.shot_id}]: Polarity array {polarity}"
                )
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, physics_method_params.times, "linear")
            dipprog_dt = interp1(
                t_ip_prog, dipprog_dt, physics_method_params.times, "linear"
            )
        except Exception as e:
            cls.logger.info(
                f"[Shot {physics_method_params.shot_id}]:Failed to get programmed plasma current parameters"
            )
            cls.logger.debug(
                f"[Shot {physics_method_params.shot_id}]:{traceback.format_exc()}"
            )
        epsoff, t_epsoff = physics_method_params.mds_conn.get_data_with_dims(
            f"ptdata('epsoff', {physics_method_params.shot_id})", tree_name="d3d"
        )
        t_epsoff = (
            t_epsoff / 1.0e3 + 0.001
        )  # [ms] -> [s] # Avoid problem with simultaneity of epsoff being triggered exactly on the last time sample
        epsoff = interp1(t_epsoff, epsoff, physics_method_params.times, "linear")
        railed_indices = np.where(np.abs(epsoff) > 0.5)
        power_supply_railed = np.zeros(len(physics_method_params.times))
        power_supply_railed[railed_indices] = 1
        indices_flattop = np.where(
            (np.abs(dipprog_dt) <= 2.0e3)
            & (np.abs(ip_prog) > 100e3)
            & (power_supply_railed != 1)
        )
        physics_method_params.times = physics_method_params.times[indices_flattop]
        return physics_method_params

    @classmethod
    def _modify_times_rampup_and_flattop_timebase(
        cls, physics_method_params: PhysicsMethodParams, **kwargs
    ):
        raise NotImplementedError("Rampup and flattop not implemented for D3D")
