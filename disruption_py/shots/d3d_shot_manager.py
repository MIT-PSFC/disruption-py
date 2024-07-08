#!/usr/bin/env python3

import traceback

import numpy as np

from disruption_py.databases.d3d_database import D3DDatabase
from disruption_py.mdsplus_integration.mds_connection import (
    MDSConnection,
    ProcessMDSConnection,
)
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.shot_manager import ShotManager
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1
from disruption_py.utils.utils import without_duplicates


class D3DShotManager(ShotManager):

    def __init__(
        self, process_database: D3DDatabase, process_mds_conn: ProcessMDSConnection
    ):
        super().__init__(
            process_database=process_database, process_mds_conn=process_mds_conn
        )

    def shot_setup(
        self, shot_id: int, shot_settings: ShotSettings, **kwargs
    ) -> ShotProps:
        """
        Sets up the shot properties for cmod.
        """

        try:
            disruption_time = self.process_database.get_disruption_time(shot_id=shot_id)
        except Exception as e:
            disruption_time = None
            self.logger.error(
                f"Failed to retreive disruption time with error {e}. Continuing as if the shot did not disrupt."
            )

        mds_conn = self.process_mds_conn.get_shot_connection(shot_id=shot_id)

        efit_tree_name = self.process_database.get_efit_tree(shot_id)

        mds_conn.add_tree_nickname_funcs(
            tree_nickname_funcs={"_efit_tree": lambda: efit_tree_name}
        )

        try:
            shot_props = self.setup_shot_props(
                shot_id=shot_id,
                mds_conn=mds_conn,
                database=self.process_database,
                disruption_time=disruption_time,
                shot_settings=shot_settings,
                tokamak=Tokamak.D3D,
                **kwargs,
            )
            return shot_props
        except Exception as e:
            self.logger.info(
                f"[Shot {shot_id}]: Caught failed to setup shot {shot_id}, cleaning up tree manager."
            )
            mds_conn.close_all_trees()
            raise e

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
        try:
            (
                ip_prog,
                t_ip_prog,
            ) = shot_props.mds_conn.get_data_with_dims(
                f"ptdata('iptipp', {shot_props.shot_id})", tree_name="d3d"
            )
            t_ip_prog = t_ip_prog / 1.0e3  # [ms] -> [s]
            polarity = np.unique(
                shot_props.mds_conn.get_data(
                    f"ptdata('iptdirect', {shot_props.shot_id})", tree_name="d3d"
                )
            )
            if len(polarity) > 1:
                cls.logger.info(
                    f"[Shot {shot_props.shot_id}]:Polarity of Ip target is not constant. Using value at first timestep."
                )
                cls.logger.debug(
                    f"[Shot {shot_props.shot_id}]: Polarity array {polarity}"
                )
                polarity = polarity[0]
            ip_prog = ip_prog * polarity
            dipprog_dt = np.gradient(ip_prog, t_ip_prog)
            ip_prog = interp1(t_ip_prog, ip_prog, shot_props.times, "linear")
            dipprog_dt = interp1(t_ip_prog, dipprog_dt, shot_props.times, "linear")
        except Exception as e:
            cls.logger.info(
                f"[Shot {shot_props.shot_id}]:Failed to get programmed plasma current parameters"
            )
            cls.logger.debug(f"[Shot {shot_props.shot_id}]:{traceback.format_exc()}")
        epsoff, t_epsoff = shot_props.mds_conn.get_data_with_dims(
            f"ptdata('epsoff', {shot_props.shot_id})", tree_name="d3d"
        )
        t_epsoff = (
            t_epsoff / 1.0e3 + 0.001
        )  # [ms] -> [s] # Avoid problem with simultaneity of epsoff being triggered exactly on the last time sample
        epsoff = interp1(t_epsoff, epsoff, shot_props.times, "linear")
        railed_indices = np.where(np.abs(epsoff) > 0.5)
        power_supply_railed = np.zeros(len(shot_props.times))
        power_supply_railed[railed_indices] = 1
        indices_flattop = np.where(
            (np.abs(dipprog_dt) <= 2.0e3)
            & (np.abs(ip_prog) > 100e3)
            & (power_supply_railed != 1)
        )
        shot_props.times = shot_props.times[indices_flattop]
        return shot_props

    @classmethod
    def _modify_times_rampup_and_flattop_timebase(cls, shot_props: ShotProps, **kwargs):
        raise NotImplementedError("Rampup and flattop not implemented for D3D")
