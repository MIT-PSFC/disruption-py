#!/usr/bin/env python3

"""
Handles output settings for retrieving and saving shot data.

This module provides classes and methods to manage various output settings
for shot data, including saving to files, databases, lists, dictionaries, and
dataframes.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Type, Union

import xarray as xr
from loguru import logger

from disruption_py.core.utils.misc import shot_log_msg
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


class DatasetOutputSetting(OutputSetting):
    """
    Outputs shot data as an xarray Dataset.
    """

    def __init__(self, filepath: str = None):
        """
        Initialize DatasetOutputSetting with an empty list of datasets.

        Parameters
        ----------
        filepath : str, optional
            The path to the file where the dataset will be saved.
        """
        self.filepath = filepath
        self.datasets = []

    def _output_shot(self, params: OutputSettingParams):
        """
        Output a single shot dataset by appending its data to the list of results.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        logger.trace(
            shot_log_msg("Appending dataset: {length} time slices"),
            shot=params.shot_id,
            length=params.result.time.shape if "time" in params.result else "no",
        )
        self.datasets.append(params.result)

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get the concatenated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        xr.Dataset
            The dataset of the concatenated results.
        """
        ds = xr.concat(self.datasets, dim="shot")
        if self.filepath and not self.filepath.endswith(".csv"):
            logger.debug(f"Saving dataset: {self.filepath}")
            ds.to_netcdf(self.filepath)
        return ds


class DataFrameOutputSetting(DatasetOutputSetting):
    """
    Output shot data as a pandas.DataFrame.
    """

    def __init__(self, *args):
        """
        Display a visible DeprecationWarning.
        """
        super().__init__(*args)
        logger.warning(
            "The pandas DataFrame output setting has been deprecated "
            "and will be removed in the next release. "
            "Please use xarray Dataset, instead."
        )

    def get_results(self, params: CompleteOutputSettingParams):
        """
        Get the concatenated results.

        Parameters
        ----------
        params : CompleteOutputSettingParams
            The parameters for output cleanup and result fetching.

        Returns
        -------
        pd.DataFrame
            The dataframe of the concatenated results.
        """
        ds = super().get_results(params)
        df = ds.to_dataframe().reset_index()
        if self.filepath:
            logger.debug(f"Saving dataframe: {self.filepath}")
            df.to_csv(self.filepath, index=False)
        return df


class SQLOutputSetting(DataFrameOutputSetting):
    """
    Stream shot data to a SQL table.
    By default, stream to the test table `disruption_warning_test`.
    """

    def __init__(
        self,
        should_update=False,
        should_override_columns: List[str] = None,
    ):
        """
        Initialize the SQLOutputSetting.

        Parameters
        ----------
        should_update : bool, optional
            Whether to update existing records (default is False).
        should_override_columns : List[str], optional
            List of columns to override in the SQL table (default is None).
        """
        super().__init__()
        self.should_update = should_update
        self.should_override_columns = should_override_columns

    def _output_shot(self, params: OutputSettingParams):
        """
        Output the current shot data to the SQL table.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters containing the result to be outputted.
        """
        super()._output_shot(params)
        df = params.result.to_dataframe().reset_index()
        logger.debug(
            shot_log_msg("Updating SQL table: {shape}"),
            shot=params.shot_id,
            shape=df.shape,
        )
        params.database.add_shot_data(
            shot_id=params.shot_id,
            shot_data=df,
            update=self.should_update,
            override_columns=self.should_override_columns,
        )


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
    ".csv": DataFrameOutputSetting,
    ".h5": DatasetOutputSetting,
    ".hdf5": DatasetOutputSetting,
    ".nc": DatasetOutputSetting,
}
# --8<-- [end:file_suffix_to_output_setting_dict]


def resolve_output_setting(output_setting: OutputSettingType) -> OutputSetting:
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
        # shortcuts
        output_setting_cls = _output_setting_mappings.get(output_setting)
        if output_setting_cls:
            return output_setting_cls()
        # extensions
        for ext, output_setting_cls in _file_suffix_to_output_setting.items():
            if output_setting.lower().endswith(ext):
                return output_setting_cls(output_setting)

    if isinstance(output_setting, list):
        return OutputSettingList(output_setting)

    raise ValueError(f"Invalid output setting: {output_setting}")
