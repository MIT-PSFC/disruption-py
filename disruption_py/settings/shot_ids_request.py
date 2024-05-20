from abc import ABC, abstractmethod
from dataclasses import dataclass
from importlib import resources
from logging import Logger
from typing import Dict, List, Type, Union

import numpy as np
import pandas as pd

import disruption_py.data
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass
class ShotIdsRequestParams:
    """Params passed by disruption_py to _get_shot_ids() method.

    Attributes
    ----------
    database : ShotDatabase
        Database object to use for getting shot ids.
        A different database connection is used by each process.
        Defaults to logbook.
    tokamak : Tokemak
        The tokamak for which the data request is made.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    database: ShotDatabase
    tokamak: Tokamak
    logger: Logger


class ShotIdsRequest(ABC):
    """ShotIdsRequest abstract class that should be inherited by all shot id request classes."""

    def get_shot_ids(self, params: ShotIdsRequestParams) -> List:
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._get_shot_ids(params)

    @abstractmethod
    def _get_shot_ids(self, params: ShotIdsRequestParams) -> List:
        """Abstract method implemented by subclasses to get shot ids for the given request params as a list.

        Parameters
        ----------
        params : ShotIdsRequestParams
            Params that can be used to determine shot ids.
        """
        pass


class IncludedShotIdsRequest(ShotIdsRequest):
    """Use the shot IDs from one of the provided data files.

    Directly passing a key from the _get_shot_ids_request_mappings dictionary as a string will
    automatically crete a new IncludedShotIdsRequest object with that data_file_name.

    Parameters
    ----------
    data_file_name : str
        The name of the datafile that should be used to retrieve shot_ids.
    """

    def __init__(self, data_file_name: str) -> List:
        with resources.path(disruption_py.data, data_file_name) as data_file:
            df = pd.read_csv(data_file, header=None)
            lst = df.values[:, 0].tolist()
            self.shot_ids = lst

    def _get_shot_ids(self, params: ShotIdsRequestParams) -> List:
        return self.shot_ids


class FileShotIdsRequest(ShotIdsRequest):
    """Use a list of shot IDs from the provided file path, this may be any file readable by pandas read_csv.

    Directly passing a file path as a string to the shot id request with the file name suffixed by txt or csv
    will automatically create a new FileShotIdsRequest object with that file path.

    Parameters
    ----------
    file_path : str
        The file path of the file that should be used for retrieving shot ids.
    column_index : int
        The index of the column that should be read. For text files, this should be 0. Defaults to 0.
    """

    def __init__(self, file_path, column_index=0):
        self.shot_ids = (
            pd.read_csv(file_path, header=None).iloc[:, column_index].values.tolist()
        )

    def _get_shot_ids(self, params: ShotIdsRequestParams) -> List:
        return self.shot_ids


class DatabaseShotIdsRequest(ShotIdsRequest):
    """Use an sql query of the database to retrieve the shot ids.

    Parameters
    ----------
    sql_query : str
        The sql query that should be used for retrieving shot ids.
    use_pandas : bool
        Whether Pandas should be used to do the sql query. Defaults to true.
    """

    def __init__(self, sql_query, use_pandas=True):
        self.sql_query = sql_query
        self.use_pandas = use_pandas

    def _get_shot_ids(self, params: ShotIdsRequestParams) -> List:
        if self.use_pandas:
            query_result_df = params.database.query(
                query=self.sql_query, use_pandas=True
            )
            return query_result_df.iloc[:, 0].tolist()
        else:
            query_result = params.database.query(query=self.sql_query, use_pandas=False)
            return [row[0] for row in query_result]


# --8<-- [start:get_shot_ids_request_dict]
_get_shot_ids_request_mappings: Dict[str, ShotIdsRequest] = {
    "d3d_paper_shotlist": IncludedShotIdsRequest("paper_shotlist.txt"),
    "d3d_train_disr": IncludedShotIdsRequest("train_disr.txt"),
    "d3d_train_nondisr": IncludedShotIdsRequest("train_nondisr.txt"),
    "cmod_test": IncludedShotIdsRequest("cmod_test.txt"),
    "cmod_non_disruptions_ids_not_blacklist": IncludedShotIdsRequest(
        "cmod_non_disruptions_ids_not_blacklist.txt"
    ),
    "cmod_non_disruptions_ids_not_blacklist_mini": IncludedShotIdsRequest(
        "cmod_non_disruptions_ids_not_blacklist_mini.txt"
    ),
}
# --8<-- [end:get_shot_ids_request_dict]

# --8<-- [start:file_suffix_to_shot_ids_request_dict]
_file_suffix_to_shot_ids_request: Dict[str, Type[ShotIdsRequest]] = {
    ".txt": FileShotIdsRequest,
    ".csv": FileShotIdsRequest,
}
# --8<-- [end:file_suffix_to_shot_ids_request_dict]

ShotIdsRequestType = Union[
    "ShotIdsRequest",
    int,
    str,
    Dict[Tokamak, "ShotIdsRequestType"],
    List["ShotIdsRequestType"],
]


def shot_ids_request_runner(shot_ids_request, params: ShotIdsRequestParams):
    if isinstance(shot_ids_request, ShotIdsRequest):
        return shot_ids_request.get_shot_ids(params)

    if isinstance(shot_ids_request, int) or (
        isinstance(shot_ids_request, str) and shot_ids_request.isdigit()
    ):
        return [shot_ids_request]

    if isinstance(shot_ids_request, np.ndarray):
        return shot_ids_request

    if isinstance(shot_ids_request, str):
        shot_ids_request_object = _get_shot_ids_request_mappings.get(
            shot_ids_request, None
        )
        if shot_ids_request_object is not None:
            return shot_ids_request_object.get_shot_ids(params)

    if isinstance(shot_ids_request, str):
        # assume that it is a file path
        for suffix, shot_ids_request_type in _file_suffix_to_shot_ids_request.items():
            if shot_ids_request.endswith(suffix):
                return shot_ids_request_type(shot_ids_request).get_shot_ids(params)

    if isinstance(shot_ids_request, dict):
        shot_ids_request = {
            map_string_to_enum(tokamak, Tokamak): shot_ids_request_mapping
            for tokamak, shot_ids_request_mapping in shot_ids_request.items()
        }
        chosen_request = shot_ids_request.get(params.tokamak, None)
        if chosen_request is not None:
            return shot_ids_request_runner(chosen_request, params)

    if isinstance(shot_ids_request, list):
        all_results = []
        for request in shot_ids_request:
            sub_result = shot_ids_request_runner(request, params)
            if sub_result is not None:
                all_results.append(sub_result)

        return [shot_id for sub_list in all_results for shot_id in sub_list]

    raise ValueError("Invalid shot id request")
