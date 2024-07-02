#!/usr/bin/env python3

import numpy as np

from disruption_py.mdsplus_integration.mds_connection import MDSConnection
from disruption_py.shots.helpers.parameter_method_params import ParameterMethodParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.parameter_methods.cmod.basic_parameter_methods import (
    BasicCmodRequests,
)
from disruption_py.shots.shot_manager import ShotManager
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.utils import without_duplicates


class CModShotManager(ShotManager):

    @classmethod
    def get_efit_tree_nickname_func(
        cls,
        shot_id: int,
        mds_conn: MDSConnection,
        disruption_time: float,
        shot_settings: ShotSettings,
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
    def _modify_times_flattop_timebase(cls, shot_props: ShotProps, **kwargs):
        parameter_method_params = ParameterMethodParams(
            mds_conn=shot_props.mds_conn,
            shot_props=shot_props,
            logger=cls.logger,
            tokamak=Tokamak.CMOD,
        )
        ip_parameters = BasicCmodRequests._get_ip_parameters(
            params=parameter_method_params
        )
        ipprog, dipprog_dt = ip_parameters["ip_prog"], ip_parameters["dipprog_dt"]
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 1e3)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.0e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            cls.logger.warning(
                f"[Shot {shot_props.shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase."
            )
            return
        shot_props.times = shot_props.times[indices_flattop]
        shot_props._cached_results.clear()  # TODO: Make this only modify the cached results for new times
        return shot_props

    @classmethod
    def _modify_times_rampup_and_flattop_timebase(cls, shot_props: ShotProps, **kwargs):
        parameter_method_params = ParameterMethodParams(
            mds_conn=shot_props.mds_conn,
            shot_props=shot_props,
            logger=cls.logger,
            tokamak=Tokamak.CMOD,
        )
        ip_parameters = BasicCmodRequests._get_ip_parameters(
            params=parameter_method_params
        )
        ipprog, dipprog_dt = ip_parameters["ip_prog"], ip_parameters["dipprog_dt"]
        indices_flattop_1 = np.where(np.abs(dipprog_dt) <= 6e4)[0]
        indices_flattop_2 = np.where(np.abs(ipprog) > 1.0e5)[0]
        indices_flattop = np.intersect1d(indices_flattop_1, indices_flattop_2)
        if len(indices_flattop) == 0:
            cls.logger.warning(
                f"[Shot {shot_props.shot_id}]:Could not find flattop timebase. Defaulting to full shot(efit) timebase."
            )
            return
        end_index = np.max(indices_flattop)
        shot_props.times = shot_props.times[:end_index]
        shot_props._cached_results.clear()  # TODO: Make this only modify the cached results for new times
        return shot_props
