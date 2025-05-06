#!/usr/bin/env python3

"""
Module for managing retrieval of shot data from a tokamak.
"""
import MDSplus.mdsExceptions
import numpy as np
from loguru import logger

from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.physics_method.runner import populate_shot
from disruption_py.core.utils.misc import shot_msg
from disruption_py.inout.mds import MDSConnection, ProcessMDSConnection
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.domain_setting import DomainSettingParams
from disruption_py.settings.nickname_setting import NicknameSettingParams
from disruption_py.settings.output_setting import OutputType
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.settings.time_setting import TimeSettingParams


class RetrievalManager:
    """
    Class for managing the retrieval of shot data from a tokamak.

    Attributes
    ----------
    tokamak : Tokamak
        The tokamak instance.
    process_database : ShotDatabase
        The SQL database
    process_mds_conn : ProcessMDSConnection
        The MDS connection
    """

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
    ) -> OutputType:
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
        pd.DataFrame, or None
            The retrieved shot data as a DataFrame, or None if an error occurred.
        """

        logger.trace(shot_msg("Starting retrieval."), shot=shot_id)

        # shot setup
        try:
            physics_method_params = self.shot_setup(
                shot_id=int(shot_id),
                retrieval_settings=retrieval_settings,
            )
        # pylint: disable-next=broad-exception-caught
        except Exception as e:
            logger.critical(shot_msg("Failed setup! {e}"), shot=shot_id, e=repr(e))
            logger.opt(exception=True).debug(shot_msg("Failed setup!"), shot=shot_id)
            return None

        # shot retrieval
        try:
            retrieved_data = populate_shot(
                retrieval_settings=retrieval_settings,
                physics_method_params=physics_method_params,
            )
        # pylint: disable-next=broad-exception-caught
        except Exception as e:
            # exceptions should be caught by runner.py
            logger.critical(shot_msg("Failed retrieval! {e}"), shot=shot_id, e=repr(e))
            logger.opt(exception=True).debug(
                shot_msg("Failed retrieval!"), shot=shot_id
            )
            if isinstance(e, MDSplus.mdsExceptions.MDSplusERROR):
                physics_method_params.mds_conn.reconnect()
            retrieved_data = None

        # shot cleanup
        try:
            self.shot_cleanup(physics_method_params)
        # pylint: disable-next=broad-exception-caught
        except Exception as e:
            logger.critical(shot_msg("Failed cleanup! {e}"), shot=shot_id, e=repr(e))
            logger.opt(exception=True).debug(shot_msg("Failed cleanup!"), shot=shot_id)
            if isinstance(e, MDSplus.mdsExceptions.MDSplusERROR):
                physics_method_params.mds_conn.reconnect()
            retrieved_data = None

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
                    )
                )
            }
        )

        physics_method_params = self.setup_physics_method_params(
            shot_id=shot_id,
            mds_conn=mds_conn,
            disruption_time=disruption_time,
            retrieval_settings=retrieval_settings,
            **kwargs,
        )
        if len(physics_method_params.times) < 2:
            raise ValueError("Pathological timebase.")
        return physics_method_params

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

        times = self._init_times(
            shot_id=shot_id,
            mds_conn=mds_conn,
            disruption_time=disruption_time,
            retrieval_settings=retrieval_settings,
        )

        physics_method_params = PhysicsMethodParams(
            shot_id=shot_id,
            tokamak=self.tokamak,
            disruption_time=disruption_time,
            mds_conn=mds_conn,
            times=times,
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
            )
        )
        if new_timebase is not None:
            physics_method_params.times = new_timebase
            # TODO: Make this only modify the cached results for new times
            physics_method_params.cached_results.clear()

        return physics_method_params

    def _init_times(
        self,
        shot_id: int,
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
            database=self.process_database,
            disruption_time=disruption_time,
            tokamak=self.tokamak,
        )
        return retrieval_settings.time_setting.get_times(setting_params)
