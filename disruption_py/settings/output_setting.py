#!/usr/bin/env python3

"""
Handles output settings for retrieving and saving shot data.

This module provides classes and methods to manage various output settings
for shot data, including saving to files, databases, lists, dictionaries, and
dataframes.
"""

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
    """
    Parameters for outputting shot results.

    Attributes
    ----------
    shot_id : int
        Shot ID.
    result : pd.DataFrame
        DataFrame of shot results.
    database : ShotDatabase
        Database connection for retrieving cache data.
    tokamak : Tokamak
        The tokamak for which results are being outputted.
    logger : Logger
        Logger instance.
    """

    shot_id: int
    result: pd.DataFrame
    database: ShotDatabase
    tokamak: Tokamak
    logger: Logger


@dataclass
class CompleteOutputSettingParams:
    """
    Parameters for output cleanup and result fetching.

    Attributes
    ----------
    tokamak : Tokamak
        The tokamak for which results are being outputted.
    logger : Logger
        Logger instance.
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
    """
    OutputSetting abstract class that should be inherited by all output setting classes.
    """

    def output_shot(self, params: OutputSettingParams):
        """
        Output a single shot based on the provided parameters.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        if hasattr(self, "tokamak_overrides"):
            if params.tokamak in self.tokamak_overrides:
                return self.tokamak_overrides[params.tokamak](params)
        return self._output_shot(params)

    @abstractmethod
    def _output_shot(self, params: OutputSettingParams):
        """
        Abstract method implemented by subclasses to handle data output for a
        single shot.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Empty method optionally overridden by subclasses to handle cleanup after
        all shots have been output. This may include closing files or other cleanup.
        Subclasses should implement this method so multiple output types can be
        used for the same data without appending to the other's outputted dataframe.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """

    @abstractmethod
    def get_results(self, params: CompleteOutputSettingParams) -> Any:
        """
        Return final output after all shots are processed.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        Any
            The final output results.
        """


class OutputSettingList(OutputSetting):
    """
    Handles a list of output settings.
    """

    def __init__(self, output_setting_list: List[OutputSettingType]):
        """
        Initialize OutputSettingList with a list of output settings.

        Parameters
        ----------
        output_setting_list : List[OutputSettingType]
            A list of output settings to handle.
        """
        self.output_setting_list = [
            resolve_output_setting(individual_setting)
            for individual_setting in output_setting_list
        ]

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot for each output setting in the list.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        return [s.output_shot(params) for s in self.output_setting_list]

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup for each output setting in the list.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        for individual_setting in self.output_setting_list:
            individual_setting.stream_output_cleanup(params)

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get results from each output setting in the list.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        List[Any]
            A list of results from each output setting.
        """
        return [s.get_results(params) for s in self.output_setting_list]


class OutputSettingDict(OutputSetting):
    """
    Handles output settings based on a dictionary of tokamaks.
    """

    def __init__(self, output_setting_dict: Dict[Tokamak, OutputSettingType]):
        """
        Initialize OutputSettingDict with a dictionary of output settings.

        Parameters
        ----------
        output_setting_dict : Dict[Tokamak, OutputSettingType]
            A dictionary mapping tokamak instances to their respective output settings.
        """
        self.output_setting_dict = {
            map_string_to_enum(tokamak, Tokamak): resolve_output_setting(setting)
            for tokamak, setting in output_setting_dict.items()
        }

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot for the specified tokamak.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        chosen_setting = self.output_setting_dict.get(params.tokamak)
        if chosen_setting is not None:
            return chosen_setting.output_shot(params)
        params.logger.warning("No output setting for tokamak %s", params.tokamak)
        return None

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup for the output setting of the specified tokamak.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        chosen_setting = self.output_setting_dict.get(params.tokamak)
        if chosen_setting:
            return chosen_setting.stream_output_cleanup(params)
        params.logger.warning("No output setting for tokamak %s", params.tokamak)
        return None

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get results from the output setting of the specified tokamak.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        Any
            The results from the output setting for the specified tokamak.
        """
        chosen_setting = self.output_setting_dict.get(params.tokamak, None)
        if chosen_setting is not None:
            return chosen_setting.get_results(params)
        params.logger.warning("No output setting for tokamak %s", params.tokamak)
        return None


class ListOutputSetting(OutputSetting):
    """
    Outputs shot data as a list of DataFrames.
    """

    def __init__(self):
        """Initialize ListOutputSetting with an empty results list."""
        self.results = []

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot by appending the result to the list.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        self.results.append(params.result)

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get the accumulated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        List[pd.DataFrame]
            The list of results.
        """
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup the results list.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        self.results = []


class DictOutputSetting(OutputSetting):
    """
    Outputs shot data as a dict of DataFrames keyed by shot number.
    """

    def __init__(self):
        """Initialize DictOutputSetting with an empty results dictionary."""
        self.results = {}

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot by storing the result in the dictionary.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        self.results[params.shot_id] = params.result

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get the accumulated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        Dict[int, pd.DataFrame]
            The dictionary of results keyed by shot number.
        """
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup the results dictionary.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        self.results = {}


class DataFrameOutputSetting(OutputSetting):
    """
    Outputs all shot data as a single DataFrame.
    """

    def __init__(self):
        """Initialize DataFrameOutputSetting with an empty DataFrame."""
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot by concatenating the result to the DataFrame.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        self.results = safe_df_concat(self.results, [params.result])

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get the accumulated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        pd.DataFrame
            The combined DataFrame of results.
        """
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup the results DataFrame.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        self.results = pd.DataFrame()


class HDF5OutputSetting(OutputSetting):
    """
    Stream outputted data to an HDF5 file. Data for each shot is stored in a table
    under the key `df_SHOTID`.
    """

    def __init__(self, filepath: str, only_output_numeric: bool = False):
        """
        Initialize HDF5OutputSetting with a file path and numeric output option.

        Parameters
        ----------
        filepath : str
            The path to the HDF5 file.
        only_output_numeric : bool, optional
            If True, only numeric data will be outputted (default is False) and
            non-numeric quantities like commit hash will be excluded.
        """
        self.filepath = filepath
        self.output_shot_count = 0
        self.only_output_numeric = only_output_numeric
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot to the HDF5 file.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
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
        """
        Get the accumulated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        pd.DataFrame
            The combined DataFrame of results.
        """
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup the results DataFrame.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        self.results = pd.DataFrame()


class CSVOutputSetting(OutputSetting):
    """
    Outputs shot data to a single CSV file.
    Not recommended when retrieving a large number of shots.
    """

    def __init__(
        self, filepath: str, flexible_columns: bool = True, clear_file: bool = True
    ):
        """
        Initialize CSVOutputSetting with a file path and options.

        Parameters
        ----------
        filepath : str
            The path to the CSV file.
        flexible_columns : bool, optional
            If True, allows for flexible columns in the CSV (default is True).
        clear_file : bool, optional
            If True, clears the file if it exists (default is True).
        """
        self.filepath = filepath
        self.flexible_columns = flexible_columns
        self.output_shot_count = 0
        if clear_file and os.path.exists(filepath):
            os.remove(filepath)
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot to the CSV file.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
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
        """
        Get the accumulated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        pd.DataFrame
            The combined DataFrame of results.
        """
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Cleanup the results DataFrame.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.
        """
        self.results = pd.DataFrame()


class BatchedCSVOutputSetting(OutputSetting):
    """
    Stream outputted data to a single CSV file in batches.
    """

    def __init__(self, filepath, batch_size=100, clear_file=True):
        """
        Initialize the BatchedCSVOutputSetting.

        Parameters
        ----------
        filepath : str
            The path to the CSV file where data will be written.
        batch_size : int, optional
            The number of records to write to the CSV file in one batch (default is 100).
        clear_file : bool, optional
            Whether to clear the file at the beginning (default is True).
        """
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
        """
        Append the current result to the batch data list and write to CSV if
        batch size is reached.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters containing the result to be outputted.
        """
        # Append the current result to the batch data list
        self.batch_data.append(params.result)

        # Check if the batch size has been reached
        if len(self.batch_data) >= self.batch_size:
            self._write_batch_to_csv()

        self.output_shot_count += 1
        self.results = safe_df_concat(self.results, [params.result])

    def _write_batch_to_csv(self):
        """
        Write the current batch of data to the CSV file.
        """
        file_exists = os.path.isfile(self.filepath)
        combined_df = safe_df_concat(pd.DataFrame(), self.batch_data)
        combined_df.to_csv(
            self.filepath, mode="a", index=False, header=(not file_exists)
        )
        self.batch_data.clear()

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Write any remaining batched data to the CSV file before returning results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for retrieving results.

        Returns
        -------
        pd.DataFrame
            The DataFrame containing the results.
        """
        # Write any remaining batched data to the CSV file before returning results
        if self.batch_data:
            self._write_batch_to_csv()
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Clean up the output stream by resetting results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for cleaning up the output stream.
        """
        self.results = pd.DataFrame()


class SQLOutputSetting(OutputSetting):
    """
    Stream outputted data to a SQL table. By default, stream to the test table:
    disruption_warning_test.
    """

    def __init__(
        self,
        should_update=False,
        should_override_columns: List[str] = None,
        table_name="disruption_warning_test",
    ):
        """
        Initialize the SQLOutputSetting.

        Parameters
        ----------
        should_update : bool, optional
            Whether to update existing records (default is False).
        should_override_columns : List[str], optional
            List of columns to override in the SQL table (default is None).
        table_name : str, optional
            The name of the SQL table to stream data to (default is
            "disruption_warning_test").
        """
        self.should_update = should_update
        self.should_override_columns = should_override_columns
        self.table_name = table_name
        self.results: pd.DataFrame = pd.DataFrame()

    def _output_shot(self, params: OutputSettingParams):
        """
        Output the current shot data to the SQL table.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters containing the result to be outputted.
        """
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
        """
        Retrieve the results stored in the SQL output setting.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for retrieving results.

        Returns
        -------
        pd.DataFrame
            The DataFrame containing the results.
        """
        return self.results

    def stream_output_cleanup(self, params: CompleteOutputSettingParams):
        """
        Clean up the output stream by resetting results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for cleaning up the output stream.
        """
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
    """
    Resolve the output setting to an OutputSetting instance.

    Parameters
    ----------
    output_setting : OutputSettingType
        The output setting to resolve, which can be an instance of OutputSetting,
        a string, a dictionary, or a list.

    Returns
    -------
    OutputSetting
        The resolved OutputSetting instance.
    """
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

    raise ValueError(f"Invalid output processor {output_setting}")
