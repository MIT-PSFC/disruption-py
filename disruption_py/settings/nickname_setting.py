#!/usr/bin/env python3

from dataclasses import dataclass
from logging import Logger

from disruption_py.core.utils.misc import without_duplicates
from disruption_py.io.mds import MDSConnection
from disruption_py.io.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class NicknameSettingParams:
    """Params passed by disruption_py to nickname trees.

    Attributes
    ----------
    sho_id : int
        the shot id for which nicknames are set
    mds_conn: MDSConnection
        MDSConnection object to access MDSPlus data.
    database : ShotDatabase
        Database connection object for tokamak sql database.
        A different database connection is used by each thread/process.
    disruption_time : float
        The time of the disruption in seconds.
    efit_tree_name : str
        The name of the efit tree requested by the user.
    tokamak : Tokamak
        The tokamak for which results are being output.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    shot_id: int
    mds_conn: MDSConnection
    database: ShotDatabase
    disruption_time: float
    efit_tree_name: str
    tokamak: Tokamak
    logger: Logger


class NicknameSetting:
    """
    A setting for getting tree nicknames.
    """

    NICKNAME_FUNC_RETRIEVERS = {
        Tokamak.D3D: "_d3d_nickname_funcs",
        Tokamak.CMOD: "_cmod_nickname_funcs",
    }

    @classmethod
    def get_nickname_funcs(cls, params: NicknameSettingParams):
        if params.tokamak in cls.NICKNAME_FUNC_RETRIEVERS:
            func_name = cls.NICKNAME_FUNC_RETRIEVERS[params.tokamak]
            return getattr(cls, func_name)(params)

        raise ValueError(
            f"No nickname function retriever for tokamak {params.tokamak}."
        )

    @classmethod
    def _cmod_nickname_funcs(cls, params: NicknameSettingParams):
        def efit_tree_nickname_func():
            efit_names_to_test = without_duplicates(
                [
                    params.efit_tree_name,
                    "analysis",
                    *[f"efit0{i}" for i in range(1, 10)],
                    *[f"efit{i}" for i in range(10, 19)],
                ]
            )

            if "efit18" in efit_names_to_test and params.disruption_time is None:
                efit_names_to_test.remove("efit18")

            for efit_name in efit_names_to_test:
                try:
                    params.mds_conn.open_tree(efit_name)
                    return efit_name
                except Exception as e:
                    cls.logger.info(
                        f"[Shot {params.shot_id}]: Failed to open efit tree {efit_name} with error {e}."
                    )
                    continue

            raise Exception(
                f"Failed to find efit tree with name {params.efit_tree_name} in shot {params.shot_id}."
            )

        return {"_efit_tree": efit_tree_nickname_func}

    @classmethod
    def _d3d_nickname_funcs(cls, params: NicknameSettingParams):
        def efit_tree_nickname_func():
            if params.efit_tree_name != "analysis":
                return params.efit_tree_name

            efit_trees = params.database.query(
                "select tree from code_rundb.dbo.plasmas where "
                f"shot = {params.shot_id} and runtag = 'DIS' and deleted = 0 order by idx",
                use_pandas=False,
            )
            if len(efit_trees) == 0:
                efit_trees = [("EFIT01",)]
            efit_tree = efit_trees[-1][0]
            return efit_tree

        return {"_efit_tree": efit_tree_nickname_func}
