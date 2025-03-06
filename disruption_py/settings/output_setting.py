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
import xarray as xr
from loguru import logger

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
    result : xr.Dataset
        Dataset of shot results.
    database : ShotDatabase
        Database connection for retrieving cache data.
    tokamak : Tokamak
        The tokamak for which results are being outputted.
    """

    shot_id: int
    result: xr.Dataset
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


def dataset_to_dataframe(ds: xr.Dataset) -> pd.DataFrame:
    """
    Convert an xarray Dataset to a pandas DataFrame.

    Parameters
    ----------
    ds : xr.Dataset
        The Dataset to convert.

    Returns
    -------
    pd.DataFrame
        The converted DataFrame.
    """
    df = ds.to_dataframe().reset_index()
    return df


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
        logger.warning(
            "No output setting for tokamak {tokamak}", tokamak=params.tokamak
        )
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
        logger.warning(
            "No output setting for tokamak {tokamak}", tokamak=params.tokamak
        )
        return None


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

        self.results = safe_df_concat(
            self.results, [dataset_to_dataframe(params.result)]
        )

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


class DatasetOutputSetting(OutputSetting):
    """
    Outputs shot data as an xarray Dataset and (optionally) save the Dataset to
    a file.
    """

    def __init__(self, filepath: str = None):
        """
        Initialize DatasetOutputSetting with an empty list of Datasets.

        Parameters
        ----------
        filepath : str, optional
            The path to the file where the dataset will be saved (default is None
            to only keep the data in memory). Accepted file extensions are .h5,
            .hdf5, and .nc.
        """
        self.filepath = filepath
        self.datasets = []

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot by concatenating the data in the Dataset.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        self.datasets.append(params.result)

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get the accumulated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        Dict[int, xr.Dataset]
            The dataset of results with dims shot and time.
        """
        ds = xr.concat(self.datasets, dim="shot")
        if self.filepath:
            ds.to_netcdf(self.filepath)
        return ds


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
        combined_df = dataset_to_dataframe(params.result)
        if self.flexible_columns:
            if file_exists:
                existing_df = pd.read_csv(self.filepath)
                combined_df = safe_df_concat(existing_df, [combined_df])

            combined_df.to_csv(self.filepath, index=False)
        else:
            params.result.to_csv(
                self.filepath, mode="a", index=False, header=(not file_exists)
            )
        self.output_shot_count += 1
        self.results = safe_df_concat(self.results, [combined_df])

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
        df = dataset_to_dataframe(params.result)
        self.batch_data.append(df)

        # Check if the batch size has been reached
        if len(self.batch_data) >= self.batch_size:
            self._write_batch_to_csv()

        self.output_shot_count += 1
        self.results = safe_df_concat(self.results, [df])

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
        df = dataset_to_dataframe(params.result)
        if not df.empty and ("shot" in df.columns):
            shot_id = df["shot"].iloc[0]
            params.database.add_shot_data(
                shot_id=shot_id,
                shot_data=df,
                update=self.should_update,
                override_columns=self.should_override_columns,
            )
        else:
            logger.warning("No shot id found in result DataFrame")
        self.results = safe_df_concat(self.results, [df])

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


# --8<-- [start:output_setting_dict]
_output_setting_mappings: Dict[str, Type[OutputSetting]] = {
    "dataframe": DataFrameOutputSetting,
    "dataset": DatasetOutputSetting,
    "pandas": DataFrameOutputSetting,
    "xarray": DatasetOutputSetting,
}
# --8<-- [end:output_setting_dict]

# --8<-- [start:file_suffix_to_output_setting_dict]
_file_suffix_to_output_setting: Dict[str, Type[OutputSetting]] = {
    ".cdf": DatasetOutputSetting,
    ".csv": BatchedCSVOutputSetting,
    ".h5": DatasetOutputSetting,
    ".hdf5": DatasetOutputSetting,
    ".nc": DatasetOutputSetting,
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
        output_setting_cls = _output_setting_mappings.get(output_setting, None)
        if output_setting_cls is not None:
            return output_setting_cls()

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
