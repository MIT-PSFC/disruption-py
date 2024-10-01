#!/usr/bin/env python3

"""
Module for managing retrieval of shot data from a tokamak.
"""

import logging

import numpy as np
import pandas as pd

from disruption_py.config import config
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.physics_method.runner import populate_shot
from disruption_py.core.utils.math import interp1
from disruption_py.core.utils.misc import get_commit_hash
from disruption_py.inout.mds import MDSConnection, ProcessMDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.cache_setting import CacheSettingParams
from disruption_py.settings.domain_setting import DomainSettingParams
from disruption_py.settings.nickname_setting import NicknameSettingParams
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.settings.time_setting import TimeSettingParams


class RetrievalManager:
    """
    Class for managing the retrieval of shot data from a tokamak.

    Attributes
    ----------
    logger : logging.Logger
        Logger for the RetrievalManager.
    tokamak : Tokamak
        The tokamak instance.
    process_database : ShotDatabase
        The SQL database
    process_mds_conn : ProcessMDSConnection
        The MDS connection
    """

    logger = logging.getLogger("disruption_py")

    def __init__(
        self,
        tokamak: Tokamak,
        process_database: ShotDatabase,
        process_mds_conn: ProcessMDSConnection,
    ):
        """
        Parameters
        ----------
        tokamak : Tokamak
            The tokamak instance.
        process_database : ShotDatabase
            The SQL database.
        process_mds_conn : ProcessMDSConnection
            The MDS connection.
        """
        self.tokamak = tokamak
        self.process_database = process_database
        self.process_mds_conn = process_mds_conn

    def get_shot_data(
        self, shot_id, retrieval_settings: RetrievalSettings
    ) -> pd.DataFrame:
        """
        Get data for a single shot. May be run across different processes.

        Parameters
        ----------
        shot_id : int
            The ID of the shot to retrieve data for.
        retrieval_settings : RetrievalSettings
            The settings for data retrieval.

        Returns
        -------
        pd.DataFrame
            The retrieved shot data as a DataFrame, or None if an error occurred.
        """
        self.logger.info("starting %s", shot_id)
        physics_method_params = self.shot_setup(
            shot_id=int(shot_id),
            retrieval_settings=retrieval_settings,
        )
        retrieved_data = populate_shot(
            retrieval_settings=retrieval_settings,
            physics_method_params=physics_method_params,
        )
        self.shot_cleanup(physics_method_params)
        self.logger.info("completed %s", shot_id)
        return retrieved_data

    def shot_setup(
        self, shot_id: int, retrieval_settings: RetrievalSettings, **kwargs
    ) -> PhysicsMethodParams:
        """
        Sets up the shot properties for the tokamak.

        Parameters
        ----------
        shot_id : int
            The ID of the shot to set up.
        retrieval_settings : RetrievalSettings
            The settings for data retrieval.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        PhysicsMethodParams
            Parameters containing MDS connection and shot information
        """

        disruption_time = self.process_database.get_disruption_time(shot_id=shot_id)

        mds_conn = self.process_mds_conn.get_shot_connection(shot_id=shot_id)

        mds_conn.add_tree_nickname_funcs(
            tree_nickname_funcs={
                "_efit_tree": lambda: retrieval_settings.efit_nickname_setting.get_tree_name(
                    NicknameSettingParams(
                        shot_id=shot_id,
                        mds_conn=mds_conn,
                        database=self.process_database,
                        disruption_time=disruption_time,
                        tokamak=self.tokamak,
                        logger=self.logger,
                    )
                )
            }
        )

        try:
            physics_method_params = self.setup_physics_method_params(
                shot_id=shot_id,
                mds_conn=mds_conn,
                disruption_time=disruption_time,
                retrieval_settings=retrieval_settings,
                **kwargs,
            )
            return physics_method_params
        except Exception as e:
            self.logger.info(
                "[Shot %s]: Caught failed to setup shot, cleaning up tree manager.",
                shot_id,
            )
            mds_conn.cleanup()
            raise e

    @classmethod
    def shot_cleanup(
        cls,
        physics_method_params: PhysicsMethodParams,
    ):
        """
        Clean up the physics method parameters.

        Parameters
        ----------
        cls : type
            The class type.
        physics_method_params : PhysicsMethodParams
            Parameters containing MDS connection and shot information.
        """
        physics_method_params.cleanup()

    def setup_physics_method_params(
        self,
        shot_id: int,
        mds_conn: MDSConnection,
        disruption_time: float,
        retrieval_settings: RetrievalSettings,
        **kwargs,
    ) -> PhysicsMethodParams:
        """
        Set up the physics method parameters for the shot.

        Parameters
        ----------
        shot_id : int
            The ID of the shot.
        mds_conn : MDSConnection
            The MDS connection for the shot.
        disruption_time : float
            The disruption time of the shot.
        retrieval_settings : RetrievalSettings
            The settings for data retrieval.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        PhysicsMethodParams
            The configured physics method parameters.
        """
        cache_data = self._retrieve_cache_data(
            shot_id=shot_id,
            retrieval_settings=retrieval_settings,
        )

        interpolation_method = interp1  # TODO: fix

        times = self._init_times(
            shot_id=shot_id,
            cache_data=cache_data,
            mds_conn=mds_conn,
            disruption_time=disruption_time,
            retrieval_settings=retrieval_settings,
        )

        pre_filled_shot_data = self._pre_fill_shot_data(
            times=times,
            cache_data=cache_data,
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
            cache_data=cache_data,
            pre_filled_shot_data=pre_filled_shot_data,
            interpolation_method=interpolation_method,
            metadata=metadata,
        )

        # Modify already existing shot properties, such as modifying timebase
        physics_method_params = self._modify_method_params_for_settings(
            physics_method_params, retrieval_settings, **kwargs
        )

        return physics_method_params

    def _modify_method_params_for_settings(
        self,
        physics_method_params: PhysicsMethodParams,
        retrieval_settings: RetrievalSettings,
        **_kwargs,
    ) -> PhysicsMethodParams:
        """
        Modify the physics method parameters based on retrieval settings.

        Parameters
        ----------
        physics_method_params : PhysicsMethodParams
            The parameters for the physics method to modify.
        retrieval_settings : RetrievalSettings
            The settings for data retrieval.
        **_kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        PhysicsMethodParams
            The modified physics method parameters.
        """
        new_timebase = retrieval_settings.domain_setting.get_domain(
            DomainSettingParams(
                physics_method_params=physics_method_params,
                tokamak=self.tokamak,
                logger=self.logger,
            )
        )
        if new_timebase is not None:
            physics_method_params.times = new_timebase
            # TODO: Make this only modify the cached results for new times
            physics_method_params.cached_results.clear()

        return physics_method_params

    def _retrieve_cache_data(
        self,
        shot_id: int,
        retrieval_settings: RetrievalSettings,
    ) -> pd.DataFrame:
        """
        Retrieve cached data for the specified shot.

        Parameters
        ----------
        shot_id : int
            The ID of the shot.
        retrieval_settings : RetrievalSettings
            The settings for data retrieval.

        Returns
        -------
        pd.DataFrame
            The cached data for the shot, or None if no cache is available.
        """
        if retrieval_settings.cache_setting is not None:
            cache_setting_params = CacheSettingParams(
                shot_id=shot_id,
                database=self.process_database,
                tokamak=self.tokamak,
                logger=self.logger,
            )
            cache_data = retrieval_settings.cache_setting.get_cache_data(
                cache_setting_params
            )
            cache_data["shot"] = cache_data["shot"].astype(int)
            cache_data = cache_data[cache_data["shot"] == shot_id]
        else:
            cache_data = None
        return cache_data

    def _init_times(
        self,
        shot_id: int,
        cache_data: pd.DataFrame,
        mds_conn: MDSConnection,
        disruption_time: float,
        retrieval_settings: RetrievalSettings,
    ) -> np.ndarray:
        """
        Initialize the timebase of the shot.

        Parameters
        ----------
        shot_id : int
            The ID of the shot.
        cache_data : pd.DataFrame
            The cached data for the shot.
        mds_conn : MDSConnection
            The MDS connection for the shot.
        disruption_time : float
            The disruption time of the shot.
        retrieval_settings : RetrievalSettings
            The settings for data retrieval.

        Returns
        -------
        np.ndarray
            The initialized timebase as a NumPy array.
        """
        setting_params = TimeSettingParams(
            shot_id=shot_id,
            mds_conn=mds_conn,
            cache_data=cache_data,
            database=self.process_database,
            disruption_time=disruption_time,
            tokamak=self.tokamak,
            logger=self.logger,
        )
        return retrieval_settings.time_setting.get_times(setting_params)

    @classmethod
    def _pre_fill_shot_data(
        cls, times: np.ndarray, cache_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Initialize the shot with data, if cached data matches the shot timebase.

        Parameters
        ----------
        cls : type
            The class type.
        times : np.ndarray
            The timebase for the shot.
        cache_data : pd.DataFrame
            The cached data for the shot.

        Returns
        -------
        pd.DataFrame
            The pre-filled shot data as a DataFrame, or None if no match is found.
        """
        if cache_data is not None:
            time_df = pd.DataFrame(times, columns=["time"])
            flagged_cache_data = cache_data.assign(merge_success_flag=1)
            timed_cache_data = pd.merge_asof(
                time_df,
                flagged_cache_data,
                on="time",
                direction="nearest",
                tolerance=config().TIME_CONST,
            )
            if not timed_cache_data["merge_success_flag"].isna().any():
                return timed_cache_data.drop(columns=["merge_success_flag"])
            return None
        return None
