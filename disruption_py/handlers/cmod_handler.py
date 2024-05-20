from typing import Callable
from disruption_py.handlers.handler import Handler
from disruption_py.shots.cmod_shot_manager import CModShotManager
from disruption_py.databases import CModDatabase
from disruption_py.utils.mappings.tokamak import Tokamak


class CModHandler(Handler):
    """Class used to retrieve MDSplus and sql data from Alcator C-Mod..

    Parameters
    ----------
    database_initializer : Callable[..., CModDatabase]
        When run returns a new database object for the handler. The function must create a new database
        connection instead of reusing an existing one, as the handler may initalize multiple connections
        across different processes. Defaults to CModDatabase.default.
    mds_connection_str : str
        The string used to connect to MDSplus using the thin client. Defaults to alcdata-new.

    Attributes
    ----------
    logger : Logger
        The logger used for disruption_py.
    database : CModDatabase
        Reference to the sql shot logbook for CMod.
    mds_connection : ProcessMDSConnection
        Reference to the MDSplus connection.

    Methods
    -------
    _get_shot_data(shot_id, sql_database=None, shot_settings)
        Static method used to get data for a single shot from CMod. Used for running across different processes.
    get_shots_data(shot_ids_request, shot_settings, num_processes)
        Instance method used to get shot data for all shots from shot_ids_request from CMod.
    """

    def __init__(
        self,
        database_initializer: Callable[..., CModDatabase] = None,
        mds_connection_str=None,
        **kwargs
    ):
        super().__init__(
            database_initializer=database_initializer or CModDatabase.default,
            mds_connection_str=mds_connection_str or "alcdata-archives",
            **kwargs
        )

    def get_shot_manager_cls(self):
        return CModShotManager

    def get_tokamak(self):
        return Tokamak.CMOD
