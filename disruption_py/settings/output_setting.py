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
from typing import Any, Dict, List, Type, Union

import pandas as pd

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
    """

    shot_id: int
    result: pd.DataFrame
    database: ShotDatabase
    tokamak: Tokamak


@dataclass
class CompleteOutputSettingParams:
    """
    Parameters for output cleanup and result fetching.

    Attributes
    ----------
    tokamak : Tokamak
        The tokamak for which results are being outputted.
    """

    tokamak: Tokamak


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
        self.columns = None

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
        # Enforce the to-be-saved combined_df to have the same column order as the first shot
        if self.columns is None:
            self.columns = combined_df.columns
        combined_df = combined_df[self.columns]
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


# --8<-- [start:output_setting_dict]
_output_setting_mappings: Dict[str, OutputSetting] = {
    "dataframe": DataFrameOutputSetting(),
    "dict": DictOutputSetting(),
}
# --8<-- [end:output_setting_dict]

# --8<-- [start:file_suffix_to_output_setting_dict]
_file_suffix_to_output_setting: Dict[str, Type[OutputSetting]] = {
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
        output_setting_object = _output_setting_mappings.get(output_setting)
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

    if isinstance(output_setting, list):
        return OutputSettingList(output_setting)

    raise ValueError(f"Invalid output processor {output_setting}")
