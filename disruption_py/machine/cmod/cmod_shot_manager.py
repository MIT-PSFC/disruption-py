#!/usr/bin/env python3

import numpy as np

from disruption_py.io.mds import MDSConnection
from disruption_py.io.sql import ShotDatabase
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.machine.cmod.basic import (
    BasicCmodRequests,
)
from disruption_py.machine.shot_manager import ShotManager
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak
from disruption_py.core.utils.misc import without_duplicates


class CModShotManager(ShotManager):

    @classmethod
    def _modify_times_flattop_timebase(
        cls, physics_method_params: PhysicsMethodParams, **kwargs
    ):
        ip_parameters = BasicCmodRequests._get_ip_parameters(
            params=physics_method_params
        )
        ipprog, dipprog_dt = ip_parameters["ip_prog"], ip_parameters["dipprog_dt"]
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 1e3)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.0e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            cls.logger.warning(
                f"[Shot {physics_method_params.shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase."
            )
            return
        physics_method_params.times = physics_method_params.times[indices_flattop]
        physics_method_params._cached_results.clear()  # TODO: Make this only modify the cached results for new times
        return physics_method_params

    @classmethod
    def _modify_times_rampup_and_flattop_timebase(
        cls, physics_method_params: PhysicsMethodParams, **kwargs
    ):
        ip_parameters = BasicCmodRequests._get_ip_parameters(
            params=physics_method_params
        )
        ipprog, dipprog_dt = ip_parameters["ip_prog"], ip_parameters["dipprog_dt"]
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 6e4)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.0e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            cls.logger.warning(
                f"[Shot {physics_method_params.shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase."
            )
            return
        end_index = np.max(indices_flattop)
        physics_method_params.times = physics_method_params.times[:end_index]
        physics_method_params._cached_results.clear()  # TODO: Make this only modify the cached results for new times
        return physics_method_params
