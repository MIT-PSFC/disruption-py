#!/usr/bin/env python3

import numpy as np

from disruption_py.io.mds import MDSConnection
from disruption_py.settings.settings import Settings
from disruption_py.machine.cmod.basic import (
    BasicCmodRequests,
)
from disruption_py.machine.shot_manager import ShotManager
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak
from disruption_py.core.utils.misc import without_duplicates


class CModShotManager(ShotManager):

    @classmethod
    def get_efit_tree_nickname_func(
        cls,
        shot_id: int,
        mds_conn: MDSConnection,
        disruption_time: float,
        shot_settings: Settings,
    ) -> None:
        def efit_tree_nickname_func():
            efit_names_to_test = without_duplicates(
                [
                    shot_settings.efit_tree_name,
                    "analysis",
                    *[f"efit0{i}" for i in range(1, 10)],
                    *[f"efit{i}" for i in range(10, 19)],
                ]
            )

            if "efit18" in efit_names_to_test and disruption_time is None:
                efit_names_to_test.remove("efit18")

            for efit_name in efit_names_to_test:
                try:
                    mds_conn.open_tree(efit_name)
                    return efit_name
                except Exception as e:
                    cls.logger.info(
                        f"[Shot {shot_id}]: Failed to open efit tree {efit_name} with error {e}."
                    )
                    continue

            raise Exception(
                f"Failed to find efit tree with name {shot_settings.efit_tree_name} in shot {shot_id}."
            )

        return efit_tree_nickname_func

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
