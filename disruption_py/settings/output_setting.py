#!/usr/bin/env python3

import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from logging import Logger
from typing import Any, Dict, List, Type, Union

import numpy as np
import pandas as pd

from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.core.utils.misc import safe_df_concat
from disruption_py.inout.sql import ShotDatabase
from disruption_py.machine.tokamak import Tokamak


@dataclass
class OutputSettingParams:
    """Params passed by disruption_py to _output_shot() method.

    Attributes
    ----------
    result : pd.DataFrame
        The DataFrame of results for a single shot.
    database : ShotDatabase
        Database object to use for getting cache data.
        A different database connection is used by each thread/process.
    tokamak : Tokamak
        The tokamak for which results are being output.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    shot_id: int
    result: pd.DataFrame
    database: ShotDatabase
    tokamak: Tokamak
    logger: Logger


@dataclass
class CompleteOutputSettingParams:
    """Params passed by disruption_py to stream_output_cleanup() and get_results methods.

    Attributes
    ----------
    tokamak : Tokamak
        The tokamak for which results are being output.
    logger : Logger
        Logger object from disruption_py to use for logging.
    """

    tokamak: Tokamak
    logger: Logger


OutputSettingType = Union[
    "OutputSetting",
    str,
    Dict[str, "OutputSettingType"],
    List["OutputSettingType"],
]


class OutputSetting(ABC):
    """OutputSetting abstract class that should be inherited by all output setting classes."""

    def output_shot(self, params: OutputSettingParams):
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._output_shot(params)

    @abstractmethod
    def _output_shot(self, params: OutputSettingParams):
        """Abstract method implemented by subclasses to handle data output for a
        single shot. This method is called by disruption_py with the shots DataFrame
        in the params object once the data has been retrieved.

        Parameters
        ----------
        params : OutputSettingParams
            Params containing the data retrieved for a shot in a DataFrame and other
            utility parameters.
        """

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """Empty method optionally overridden by subclasses to handle cleanup after
        all shots have been output. This may include closing files or other cleanup.
        Subclasses should implement this method so multiple output types can be
        used for the same data without appending to the other's outputted dataframe.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            Utility parameters such as the tokamak and logger.
        """

    @abstractmethod
    def get_results(self, params: CompleteOutputSettingParams) -> Any:
        """Abstract method implemented by subclasses to handle the output of the data from
        calls to `get_shots_data`. This method is called by disruption_py once `output_shot()` has been
        called for all shot ids.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            Utility parameters such as the tokamak and logger.

        Returns
        -------
        Any
            The desired output of the call to `get_shots_data` potentially containing the data for all shots,
            some aggregation of that data, or nothing.
        """


class OutputSettingList(OutputSetting):
    """
    Utility class that is automatically used when a list is passed as the `output_setting` parameter in `RetrievalSettings.

    All listed output types will be output to in the order listed.
    Similarly, results will be returned in the order listed.

    Parameters
    ----------
    output_setting_list : list[OutputSettingType]
        A python list of any other output setting option that can be passed as the `output_setting` parameter in `RetrievalSettings`.
        Any other option passable to the `output_setting` parameter in `RetrievalSettings` may be used.
    """

    def __init__(self, output_setting_list: List[OutputSettingType]):
        self.output_setting_list = [
            resolve_output_setting(individual_setting)
            for individual_setting in output_setting_list
        ]

    def _output_shot(self, params: OutputSettingParams):
        all_results = []
        for individual_setting in self.output_setting_list:
            sub_result = individual_setting.output_shot(params)
            all_results.append(sub_result)

        return all_results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        for individual_setting in self.output_setting_list:
            individual_setting.stream_output_cleanup(params)

    def get_results(self, params: CompleteOutputSettingParams):
        return [
            individual_setting.get_results(params)
            for individual_setting in self.output_setting_list
        ]


class OutputSettingDict(OutputSetting):
    """
    Utility class that is automatically used when a dictionary is passed as the `output_setting` parameter in `RetrievalSettings.

    Parameters
    ----------
    output_setting_dict : dict[Tokamak, OutputSettingType]
        A dictionary mapping tokamak type strings to the desired `OutputSettingType` for that tokamak.  E.g. `{'cmod': 'list'}`.
    """

    def __init__(self, output_setting_dict: Dict[Tokamak, OutputSettingType]):
        resolved_output_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_output_setting(
                individual_setting
            )
            for tokamak, individual_setting in output_setting_dict.items()
        }
        self.output_setting_dict = resolved_output_setting_dict

    def _output_shot(self, params: OutputSettingParams):
        chosen_setting = self.output_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.output_shot(params)
        params.logger.warning("No output setting for tokamak %s", params.tokamak)
        return None

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        chosen_setting = self.output_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.stream_output_cleanup(params)
        params.logger.warning("No output setting for tokamak %s", params.tokamak)
        return None

    def get_results(self, params: CompleteOutputSettingParams):
        chosen_setting = self.output_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_results(params)
        params.logger.warning("No output setting for tokamak %s", params.tokamak)
        return None


class ListOutputSetting(OutputSetting):
    """
    Output all retrieved shot data as a list of DataFrames, once retrieval completes.
    """

    def __init__(self):
        self.results = []

    def _output_shot(self, params: OutputSettingParams):
        self.results.append(params.result)

    def get_results(self, params: CompleteOutputSettingParams):
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = []


class DictOutputSetting(OutputSetting):
    """
    Output all retrieved shot data as a dict of DataFrames with the keys being shot numbers.
    """

    def __init__(self):
        self.results = {}

    def _output_shot(self, params: OutputSettingParams):
        self.results[params.shot_id] = params.result

    def get_results(self, params: CompleteOutputSettingParams):
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = {}


class DataFrameOutputSetting(OutputSetting):
    """
    Output all retrieved shot data as a list of DataFrames, once retrieval completes.
    """

    def __init__(self):
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        self.results = safe_df_concat(self.results, [params.result])

    def get_results(self, params: CompleteOutputSettingParams):
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = pd.DataFrame()


class HDF5OutputSetting(OutputSetting):
    """
    Stream outputted data to an HDF5 file. Data for each shot is stored in a table
    under the key `df_SHOTID`.
    """

    def __init__(self, filepath, only_output_numeric=False):
        """
        Params:
        filepath: str
            hdf5 file path
        only_output_numeric: bool
            (default False) whether the hdf5 file should exclude all non-numeric
            quantities (e.g. commit hash).
        """
        self.filepath = filepath
        self.output_shot_count = 0
        self.only_output_numeric = only_output_numeric
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        mode = "a" if self.output_shot_count > 0 else "w"

        if self.only_output_numeric:
            output_result = params.result.select_dtypes([np.number])
        else:
            output_result = params.result

        output_result.to_hdf(
            self.filepath,
            key=f"df_{params.shot_id}",
            format="table",
            complib="blosc",
            mode=mode,
        )
        self.output_shot_count += 1
        self.results = safe_df_concat(self.results, [params.result])

    def get_results(self, params: CompleteOutputSettingParams):
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = pd.DataFrame()


class CSVOutputSetting(OutputSetting):
    """
    Stream outputted data to a single csv file.

    Not recommended when retrieving a large number of shots.
    """

    def __init__(self, filepath, flexible_columns=True, clear_file=True):
        self.filepath = filepath
        self.flexible_columns = flexible_columns
        self.output_shot_count = 0
        if clear_file is True and os.path.exists(filepath):
            os.remove(filepath)
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        file_exists = os.path.isfile(self.filepath)
        if self.flexible_columns:
            if file_exists:
                existing_df = pd.read_csv(self.filepath)
                combined_df = safe_df_concat(existing_df, [params.result])
            else:
                combined_df = params.result

            combined_df.to_csv(self.filepath, index=False)
        else:
            params.result.to_csv(
                self.filepath, mode="a", index=False, header=(not file_exists)
            )
        self.output_shot_count += 1
        self.results = safe_df_concat(self.results, [params.result])

    def get_results(self, params: CompleteOutputSettingParams):
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = pd.DataFrame()


class BatchedCSVOutputSetting(OutputSetting):
    """
    Stream outputted data to a single csv file in batches.
    """

    def __init__(self, filepath, batch_size=100, clear_file=True):
        self.filepath = filepath
        self.batch_size = batch_size
        self.clear_file = clear_file
        self.batch_data = []  # Initialize an empty list to hold batched data
        self.output_shot_count = 0

        # Clear the file at the beginning if required
        if self.clear_file and os.path.exists(filepath):
            os.remove(filepath)

        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        # Append the current result to the batch data list
        self.batch_data.append(params.result)

        # Check if the batch size has been reached
        if len(self.batch_data) >= self.batch_size:
            self._write_batch_to_csv()

        self.output_shot_count += 1
        self.results = safe_df_concat(self.results, [params.result])

    def _write_batch_to_csv(self):
        file_exists = os.path.isfile(self.filepath)
        combined_df = safe_df_concat(pd.DataFrame(), self.batch_data)
        combined_df.to_csv(
            self.filepath, mode="a", index=False, header=(not file_exists)
        )
        self.batch_data.clear()

    def get_results(self, params: CompleteOutputSettingParams):
        # Write any remaining batched data to the CSV file before returning results
        if self.batch_data:
            self._write_batch_to_csv()
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = pd.DataFrame()


class SQLOutputSetting(OutputSetting):
    """
    Stream outputted data to disruption_warning or similar SQL table. By default,
    stream to the test table: disruption_warning_test.
    """

    def __init__(
        self,
        should_update=False,
        should_override_columns: List[str] = None,
        table_name="disruption_warning_test",
    ):
        self.should_update = should_update
        self.should_override_columns = should_override_columns
        self.table_name = table_name
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        if not params.result.empty and ("shot" in params.result.columns):
            shot_id = params.result["shot"].iloc[0]
            params.database.add_shot_data(
                shot_id=shot_id,
                shot_data=params.result,
                update=self.should_update,
                override_columns=self.should_override_columns,
            )
        else:
            params.logger.warning("No shot id found in result DataFrame")
        self.results = safe_df_concat(self.results, [params.result])

    def get_results(self, params: CompleteOutputSettingParams) -> Any:
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        self.results = pd.DataFrame()


# --8<-- [start:output_setting_dict]
_output_setting_mappings: Dict[str, OutputSetting] = {
    "list": ListOutputSetting(),
    "dataframe": DataFrameOutputSetting(),
    "dict": DictOutputSetting(),
}
# --8<-- [end:output_setting_dict]

# --8<-- [start:file_suffix_to_output_setting_dict]
_file_suffix_to_output_setting: Dict[str, Type[OutputSetting]] = {
    ".h5": HDF5OutputSetting,
    ".hdf5": HDF5OutputSetting,
    ".csv": BatchedCSVOutputSetting,
}
# --8<-- [end:file_suffix_to_output_setting_dict]


def resolve_output_setting(
    output_setting: OutputSettingType,
) -> OutputSetting:
    if isinstance(output_setting, OutputSetting):
        return output_setting

    if isinstance(output_setting, str):
        output_setting_object = _output_setting_mappings.get(output_setting, None)
        if output_setting_object is not None:
            return output_setting_object

    if isinstance(output_setting, str):
        # assume that it is a file path
        for (
            suffix,
            output_setting_type,
        ) in _file_suffix_to_output_setting.items():
            if output_setting.endswith(suffix):
                return output_setting_type(output_setting)

    if isinstance(output_setting, dict):
        return OutputSettingDict(output_setting)

    if isinstance(output_setting, list):
        return OutputSettingList(output_setting)

    raise ValueError(f"Invalid output processror {output_setting}")
