#!/usr/bin/env python3

import logging
from abc import ABC, abstractmethod
import traceback

import numpy as np
import pandas as pd

from disruption_py.databases.database import ShotDatabase
from disruption_py.mdsplus_integration.mds_connection import (
    MDSConnection,
    ProcessMDSConnection,
)
from disruption_py.settings.enum_options import SignalDomain
from disruption_py.settings.existing_data_request import ExistingDataRequestParams
from disruption_py.settings.set_times_request import SetTimesRequestParams
from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.populate_shot import populate_shot
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.command_utils import get_commit_hash
from disruption_py.utils.constants import TIME_CONST
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.math_utils import interp1


class ShotManager(ABC):
    logger = logging.getLogger("disruption_py")

    def __init__(
        self, process_database: ShotDatabase, process_mds_conn: ProcessMDSConnection
    ):
        self.process_database = process_database
        self.process_mds_conn = process_mds_conn

    @classmethod
    def get_shot_data(
	    cls, shot_id, shot_manager: 'ShotManager', shot_settings: ShotSettings
    ) -> pd.DataFrame:
        """
        Get data for a single shot. May be run across different processes.
        """
        cls.logger.info(f"starting {shot_id}")
        try:
            shot_props = shot_manager.shot_setup(
                shot_id=int(shot_id),
                shot_settings=shot_settings,
            )
            retrieved_data = shot_manager.shot_data_retrieval(
                shot_props=shot_props, shot_settings=shot_settings
            )
            shot_manager.shot_cleanup(shot_props)
            cls.logger.info(f"completed {shot_id}")
            return retrieved_data
        except Exception as e:
            cls.logger.warning(
                f"[Shot {shot_id}]: fatal error {traceback.format_exc()}"
            )
            cls.logger.error(f"failed {shot_id} with error {e}")
            return None
    
    @classmethod
    @abstractmethod
    def _modify_times_flattop_timebase(
        cls, shot_props: ShotProps, **kwargs
    ) -> ShotProps:
        pass

    @classmethod
    @abstractmethod
    def _modify_times_rampup_and_flattop_timebase(
        cls, shot_props: ShotProps, **kwargs
    ) -> ShotProps:
        pass

    @abstractmethod
    def shot_setup(
        self, shot_id: int, shot_settings: ShotSettings, **kwargs
    ) -> ShotProps:
        pass

    def shot_data_retrieval(self, shot_props: ShotProps, shot_settings: ShotSettings):
        shot_data_request_params = ShotDataRequestParams(
            mds_conn=shot_props.mds_conn,
            shot_props=shot_props,
            logger=self.logger,
            tokamak=shot_props.tokamak,
        )
        return populate_shot(
            shot_settings=shot_settings, params=shot_data_request_params
        )

    @classmethod
    def shot_cleanup(
        cls,
        shot_props: ShotProps,
    ):
        shot_props.cleanup()

    def setup_shot_props(
        self,
        shot_id: int,
        mds_conn: MDSConnection,
        disruption_time: float,
        shot_settings: ShotSettings,
        tokamak: Tokamak,
        **kwargs,
    ) -> ShotProps:

        existing_data = self._retrieve_existing_data(
            shot_id=shot_id,
            tokamak=tokamak,
            shot_settings=shot_settings,
        )

        interpolation_method = interp1  # TODO: fix

        times = self._init_times(
            shot_id=shot_id,
            existing_data=existing_data,
            mds_conn=mds_conn,
            tokamak=tokamak,
            disruption_time=disruption_time,
            shot_settings=shot_settings,
        )

        pre_filled_shot_data = self._pre_fill_shot_data(
            times=times,
            existing_data=existing_data,
        )

        metadata = {
            "labels": {},
            "commit_hash": get_commit_hash(),
            "timestep": {},
            "duration": {},
            "description": "",
            "disrupted": 100,  # TODO: Fix
        }

        shot_props = ShotProps(
            shot_id=shot_id,
            tokamak=tokamak,
            disruption_time=disruption_time,
            mds_conn=mds_conn,
            times=times,
            existing_data=existing_data,
            pre_filled_shot_data=pre_filled_shot_data,
            interpolation_method=interpolation_method,
            metadata=metadata,
        )

        # modify already existing shot props, such as modifying timebase
        shot_props = self._modify_shot_props_for_settings(
            shot_props, shot_settings, **kwargs
        )

        return shot_props

    def _modify_shot_props_for_settings(
        self, shot_props: ShotProps, shot_settings: ShotSettings, **kwargs
    ) -> ShotProps:
        if shot_settings.signal_domain is SignalDomain.FLATTOP:
            shot_props = self._modify_times_flattop_timebase(shot_props)
        elif shot_settings.signal_domain is SignalDomain.RAMP_UP_AND_FLATTOP:
            shot_props = self._modify_times_rampup_and_flattop_timebase(shot_props)
        if shot_props is None:
            raise ValueError(f"Shot_props set to None in modify_shot_props()")

        return shot_props

    def _retrieve_existing_data(
        self,
        shot_id: int,
        tokamak: Tokamak,
        shot_settings: ShotSettings,
    ) -> pd.DataFrame:
        if shot_settings.existing_data_request is not None:
            existing_data_request_params = ExistingDataRequestParams(
                shot_id=shot_id,
                database=self.process_database,
                tokamak=tokamak,
                logger=self.logger,
            )
            existing_data = shot_settings.existing_data_request.get_existing_data(
                existing_data_request_params
            )
            existing_data["shot"] = existing_data["shot"].astype(int)
            existing_data = existing_data[existing_data["shot"] == shot_id]
        else:
            existing_data = None
        return existing_data

    def _init_times(
        self,
        shot_id: int,
        existing_data: pd.DataFrame,
        mds_conn: MDSConnection,
        tokamak: Tokamak,
        disruption_time: float,
        shot_settings: ShotSettings,
    ) -> np.ndarray:
        """
        Initialize the timebase of the shot.
        """
        request_params = SetTimesRequestParams(
            shot_id=shot_id,
            mds_conn=mds_conn,
            existing_data=existing_data,
            database=self.process_database,
            disruption_time=disruption_time,
            tokamak=tokamak,
            logger=self.logger,
        )
        return shot_settings.set_times_request.get_times(request_params)

    @classmethod
    def _pre_fill_shot_data(
        cls, times: np.ndarray, existing_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Intialize the shot with data, if existing data matches the shot timebase.
        """
        if existing_data is not None:
            time_df = pd.DataFrame(times, columns=["time"])
            flagged_existing_data = existing_data.assign(merge_success_flag=1)
            timed_existing_data = pd.merge_asof(
                time_df,
                flagged_existing_data,
                on="time",
                direction="nearest",
                tolerance=TIME_CONST,
            )
            if not timed_existing_data["merge_success_flag"].isna().any():
                return timed_existing_data.drop(columns=["merge_success_flag"])
            else:
                return None
        else:
            return None
