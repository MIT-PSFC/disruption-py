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

import pandas as pd
import xarray as xr
from loguru import logger

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
    tokamak : Tokamak
        The tokamak for which results are being outputted.
    """

    shot_id: int
    result: xr.Dataset
    tokamak: Tokamak


OutputSettingType = Union["OutputSetting", str, List["OutputSettingType"]]


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
    def get_results(self) -> Any:
        """
        Return final output after all shots are processed.

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

    def get_results(self):
        """
        Get results from each output setting in the list.

        Returns
        -------
        List[Any]
            A list of results from each output setting.
        """
        return [s.get_results() for s in self.output_setting_list]


class DictOutputSetting(OutputSetting):
    """
    Outputs data as a dictionary of Datasets.
    """

    def __init__(self):
        """
        Initialize an empty dictionary.
        """
        self.results = {}

    def _output_shot(self, params: OutputSettingParams):
        """
        Store a single Dataset in the dictionary.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        self.results[params.shot_id] = params.result

    def get_results(self) -> Dict[int, xr.Dataset]:
        """
        Get the resulting dictionary.

        Returns
        -------
        Dict[int, xr.Dataset]
            The dictionary of results, with shots as keys.
        """
        return self.results


class DatasetOutputSetting(DictOutputSetting):
    """
    Outputs data as a single Dataset.
    """

    def get_results(self) -> xr.Dataset:
        """
        Get the resulting Dataset.

        Returns
        -------
        xr.Dataset
            The Dataset containing the results.
        """
        logger.debug("Concatenating {tot} results.", tot=len(self.results))
        return xr.concat(self.results.values(), dim="idx")


class DataTreeOutputSetting(DictOutputSetting):
    """
    Outputs data as a Datatree of Datasets.
    """

    def get_results(self) -> xr.DataTree:
        """
        Get the resulting DataTree.

        Returns
        -------
        xr.DataTree
            The DataTree containing the results, with shots as keys.
        """
        logger.debug("Branching {tot} results.", tot=len(self.results))
        return xr.DataTree.from_dict(self.results)


class DataFrameOutputSetting(DatasetOutputSetting):
    """
    Outputs data as a DataFrame.
    """

    def get_results(self) -> pd.DataFrame:
        """
        Get the resulting DataFrame.

        Returns
        -------
        pd.DataFrame
            The combined DataFrame of results.
        """
        base = ["shot", "time"]
        df = super().get_results().to_dataframe()
        cols = base + [c for c in sorted(df.columns) if c not in base]
        return df[cols]


# --8<-- [start:output_setting_dict]
_output_setting_mappings: Dict[str, OutputSetting] = {
    "dataframe": DataFrameOutputSetting(),
    "dataset": DatasetOutputSetting(),
    "datatree": DataTreeOutputSetting(),
    "df": DataFrameOutputSetting(),
    "dict": DictOutputSetting(),
    "ds": DatasetOutputSetting(),
    "dt": DataTreeOutputSetting(),
    "pandas": DataFrameOutputSetting(),
    "pd": DataFrameOutputSetting(),
    "xarray": DatasetOutputSetting(),
    "xr": DatasetOutputSetting(),
}
# --8<-- [end:output_setting_dict]

# --8<-- [start:file_suffix_to_output_setting_dict]
_file_suffix_to_output_setting: Dict[str, Type[OutputSetting]] = {
    ".cdf": DatasetOutputSetting,
    ".csv": DataFrameOutputSetting,
    ".hdf5": DatasetOutputSetting,
    ".h5": DatasetOutputSetting,
    ".nc": DatasetOutputSetting,
    "/": DictOutputSetting,
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
        # check shortcuts
        output_setting_object = _output_setting_mappings.get(output_setting)
        if output_setting_object is not None:
            return output_setting_object
        # check suffixes
        for suffix, output_setting_type in _file_suffix_to_output_setting.items():
            if output_setting.endswith(suffix):
                return output_setting_type(output_setting)

    if isinstance(output_setting, list):
        return OutputSettingList(output_setting)

    raise ValueError(f"Invalid output processor {output_setting}")
