#!/usr/bin/env python3

import logging
from abc import ABC, abstractmethod
import traceback

import numpy as np
import pandas as pd

from disruption_py.io.sql import ShotDatabase
from disruption_py.io.mds import (
    MDSConnection,
    ProcessMDSConnection,
)
from disruption_py.settings.settings import SignalDomain
from disruption_py.settings.input_setting import InputSettingParams
from disruption_py.settings.time_setting import TimeSettingParams
from disruption_py.settings.settings import Settings
from disruption_py.core.physics_method.runner import populate_shot
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.misc import get_commit_hash
from disruption_py.utils.constants import TIME_CONST
from disruption_py.machine.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1


class ShotManager(ABC):
    logger = logging.getLogger("disruption_py")

    def __init__(
        self,
        tokamak: Tokamak,
        process_database: ShotDatabase,
        process_mds_conn: ProcessMDSConnection,
    ):
        self.tokamak = tokamak
        self.process_database = process_database
        self.process_mds_conn = process_mds_conn

    def get_shot_data(self, shot_id, shot_settings: Settings) -> pd.DataFrame:
        """
        Get data for a single shot. May be run across different processes.
        """
        self.logger.info(f"starting {shot_id}")
        try:
            physics_method_params = self.shot_setup(
                shot_id=int(shot_id),
                shot_settings=shot_settings,
            )
            retrieved_data = populate_shot(
                shot_settings=shot_settings,
                physics_method_params=physics_method_params,
            )
            self.shot_cleanup(physics_method_params)
            self.logger.info(f"completed {shot_id}")
            return retrieved_data
        except Exception as e:
            self.logger.warning(
                f"[Shot {shot_id}]: fatal error {traceback.format_exc()}"
            )
            self.logger.error(f"failed {shot_id} with error {e}")
            return None

    @classmethod
    @abstractmethod
    def _modify_times_flattop_timebase(
        cls, physics_method_params: PhysicsMethodParams, **kwargs
    ) -> PhysicsMethodParams:
        pass

    @classmethod
    @abstractmethod
    def _modify_times_rampup_and_flattop_timebase(
        cls, physics_method_params: PhysicsMethodParams, **kwargs
    ) -> PhysicsMethodParams:
        pass

    def shot_setup(
        self, shot_id: int, shot_settings: Settings, **kwargs
    ) -> PhysicsMethodParams:
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

        mds_conn.add_tree_nickname_funcs(
            tree_nickname_funcs={
                "_efit_tree": self.get_efit_tree_nickname_func(
                    shot_id=shot_id,
                    mds_conn=mds_conn,
                    disruption_time=disruption_time,
                    shot_settings=shot_settings,
                )
            }
        )

        try:
            physics_method_params = self.setup_physics_method_params(
                shot_id=shot_id,
                mds_conn=mds_conn,
                disruption_time=disruption_time,
                shot_settings=shot_settings,
                **kwargs,
            )
            return physics_method_params
        except Exception as e:
            self.logger.info(
                f"[Shot {shot_id}]: Caught failed to setup shot {shot_id}, cleaning up tree manager."
            )
            mds_conn.cleanup()
            raise e

    @classmethod
    def shot_cleanup(
        cls,
        physics_method_params: PhysicsMethodParams,
    ):
        physics_method_params.cleanup()

    def setup_physics_method_params(
        self,
        shot_id: int,
        mds_conn: MDSConnection,
        disruption_time: float,
        shot_settings: Settings,
        **kwargs,
    ) -> PhysicsMethodParams:

        input_data = self._retrieve_input_data(
            shot_id=shot_id,
            shot_settings=shot_settings,
        )

        interpolation_method = interp1  # TODO: fix

        times = self._init_times(
            shot_id=shot_id,
            input_data=input_data,
            mds_conn=mds_conn,
            disruption_time=disruption_time,
            shot_settings=shot_settings,
        )

        pre_filled_shot_data = self._pre_fill_shot_data(
            times=times,
            input_data=input_data,
        )

        metadata = {
            "labels": {},
            "commit_hash": get_commit_hash(),
            "timestep": {},
            "duration": {},
            "description": "",
            "disrupted": 100,  # TODO: Fix
        }

        physics_method_params = PhysicsMethodParams(
            shot_id=shot_id,
            tokamak=self.tokamak,
            disruption_time=disruption_time,
            mds_conn=mds_conn,
            times=times,
            input_data=input_data,
            pre_filled_shot_data=pre_filled_shot_data,
            interpolation_method=interpolation_method,
            metadata=metadata,
        )

        # modify already existing shot props, such as modifying timebase
        physics_method_params = self._modify_method_params_for_settings(
            physics_method_params, shot_settings, **kwargs
        )

        return physics_method_params

    def _modify_method_params_for_settings(
        self,
        physics_method_params: PhysicsMethodParams,
        shot_settings: Settings,
        **kwargs,
    ) -> PhysicsMethodParams:
        if shot_settings.signal_domain is SignalDomain.FLATTOP:
            physics_method_params = self._modify_times_flattop_timebase(
                physics_method_params
            )
        elif shot_settings.signal_domain is SignalDomain.RAMP_UP_AND_FLATTOP:
            physics_method_params = self._modify_times_rampup_and_flattop_timebase(
                physics_method_params
            )
        if physics_method_params is None:
            raise ValueError(
                f"physics_method_params set to None in _modify_method_params_for_settings()"
            )

        return physics_method_params

    def _retrieve_input_data(
        self,
        shot_id: int,
        shot_settings: Settings,
    ) -> pd.DataFrame:
        if shot_settings.input_setting is not None:
            input_setting_params = InputSettingParams(
                shot_id=shot_id,
                database=self.process_database,
                tokamak=self.tokamak,
                logger=self.logger,
            )
            input_data = shot_settings.input_setting.get_input_data(
                input_setting_params
            )
            input_data["shot"] = input_data["shot"].astype(int)
            input_data = input_data[input_data["shot"] == shot_id]
        else:
            input_data = None
        return input_data

    def _init_times(
        self,
        shot_id: int,
        input_data: pd.DataFrame,
        mds_conn: MDSConnection,
        disruption_time: float,
        shot_settings: Settings,
    ) -> np.ndarray:
        """
        Initialize the timebase of the shot.
        """
        setting_params = TimeSettingParams(
            shot_id=shot_id,
            mds_conn=mds_conn,
            input_data=input_data,
            database=self.process_database,
            disruption_time=disruption_time,
            tokamak=self.tokamak,
            logger=self.logger,
        )
        return shot_settings.time_setting.get_times(setting_params)

    @classmethod
    def _pre_fill_shot_data(
        cls, times: np.ndarray, input_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Intialize the shot with data, if input data matches the shot timebase.
        """
        if input_data is not None:
            time_df = pd.DataFrame(times, columns=["time"])
            flagged_input_data = input_data.assign(merge_success_flag=1)
            timed_input_data = pd.merge_asof(
                time_df,
                flagged_input_data,
                on="time",
                direction="nearest",
                tolerance=TIME_CONST,
            )
            if not timed_input_data["merge_success_flag"].isna().any():
                return timed_input_data.drop(columns=["merge_success_flag"])
            else:
                return None
        else:
            return None
