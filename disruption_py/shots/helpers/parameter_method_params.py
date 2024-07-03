from disruption_py.io.mds import MDSConnection
from disruption_py.shots.shot_props import ShotProps
from disruption_py.machine.tokamak import Tokamak


from dataclasses import dataclass
from logging import Logger


@dataclass
class ParameterMethodParams:
    """Params passed by disruption_py to decorated methods.

    Attributes
    ----------
    mds_conn : ShotConnection
        The shot connection object containing the connection to MDSPlus. The same as shot_props.mds_conn.
    shot_props : ShotProps
                A reference to the shot props object containing the setup information, such as the shot id,
        timebase, and disruption time, for the shot data retrieval from MDSPlus.
    tokamak : Tokemak
        The tokamak for which data is being retrieved.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    mds_conn: MDSConnection
    shot_props: ShotProps
    tokamak: Tokamak
    logger: Logger
