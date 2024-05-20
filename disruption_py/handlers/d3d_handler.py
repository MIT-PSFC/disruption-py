from typing import Callable

from disruption_py.databases import D3DDatabase
from disruption_py.handlers.handler import Handler
from disruption_py.shots.d3d_shot_manager import D3DShotManager
from disruption_py.utils.mappings.tokamak import Tokamak


class D3DHandler(Handler):
    """
    Class used to retrieve MDSplus and sql data from D3D.

    Parameters
    ----------
    database_initializer : Callable[..., D3DDatabase]
        When run returns a new database object for the handler. The function must create a new database
        connection instead of reusing an existing one, as the handler may initalize multiple connections
        across different processes. Defaults to D3DDatabase.default.
    mds_connection_str : str
        The string used to connect to MDSplus using the thin client. Defaults to alcdata-new.

    Attributes
    ----------
    logger : Logger
        The logger used for disruption_py.
    database : D3DDatabase
        Reference to the sql shot logbook for D3D.
        mds_connection : ProcessMDSConnection
        Reference to the MDSplus connection.

    Methods
    -------
    _get_shot_data(shot_id, sql_database=None, shot_settings)
        Static method used to get data for a single shot from D3D. Used for running across different processes.
    get_shots_data(shot_ids_request, shot_settings, num_processes)
        Instance method used to get shot data for all shots from shot_ids_request from D3D.
    """

    def __init__(
        self,
        database_initializer: Callable[..., D3DDatabase] = None,
        mds_connection_str=None,
        **kwargs
    ):
        super().__init__(
            database_initializer=database_initializer or D3DDatabase.default,
            mds_connection_str=mds_connection_str or "atlas",
            **kwargs
        )

    def get_shot_manager_cls(self):
        return D3DShotManager

    def get_tokamak(self):
        return Tokamak.D3D
