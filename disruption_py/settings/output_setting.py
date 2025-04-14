#!/usr/bin/env python3

"""
Handles output settings for retrieving and saving shot data.

This module provides classes and methods to manage various output settings
for shot data, including saving to files, databases, lists, dictionaries, and
dataframes.
"""

import os
import sys
import tempfile
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Type, Union

import pandas as pd
import xarray as xr
from loguru import logger

from disruption_py.core.utils.misc import get_temporary_folder
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
OutputDictType = Dict[str, xr.Dataset]
OutputSingleType = xr.Dataset | xr.DataTree | pd.DataFrame
OutputType = OutputDictType | OutputSingleType


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
                self.tokamak_overrides[params.tokamak](params)
                return
        self._output_shot(params)

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
    def get_results(self) -> OutputType:
        """
        Return final output after all shots are processed.

        Returns
        -------
        Any
            The final output results.
        """

    @abstractmethod
    def to_disk(self) -> str:
        """
        Save final output to disk.
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
        _ = [s.output_shot(params) for s in self.output_setting_list]

    def get_results(self) -> List[OutputType]:
        """
        Get results from each output setting in the list.

        Returns
        -------
        List[OutputType]
            A list of results from each output setting.
        """
        return [s.get_results() for s in self.output_setting_list]

    def to_disk(self) -> List[str]:
        """
        Save each OutputSettingList to disk.
        """
        return [s.to_disk() for s in self.output_setting_list]


class DictOutputSetting(OutputSetting):
    """
    Outputs data as a dictionary of Datasets.
    """

    def __init__(self, path: str | bool = True):
        """
        Initialize empty DictOutputSetting.

        Parameters
        ----------
        path : str | bool, default = True
            The path for writing results to disk.
            If True, a temporary location will be used.
            If False, no results are written to disk.
        """
        self.results = {}
        if path is True:
            path = os.path.join(get_temporary_folder(), "output")
            if os.path.exists(path) or "pytest" in sys.modules:
                path = tempfile.mkdtemp(dir=get_temporary_folder(), prefix="output-")
        self.path = path

    def _output_shot(self, params: OutputSettingParams):
        """
        Store a single result in the dictionary.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """
        self.results[str(params.shot_id)] = params.result

    def get_results(self) -> OutputDictType:
        """
        Get the resulting dictionary.

        Returns
        -------
        Dict[int, xr.Dataset]
            The dictionary of results, with shots as keys.
        """
        return self.results

    def to_disk(self) -> str:
        """
        Save all resulting Datasets into a folder.
        """

        if not self.path:
            return ""
        if os.path.exists(self.path):
            if not os.path.isdir(self.path):
                raise FileExistsError(f"Path already exists! {self.path}")
            if os.listdir(self.path):
                logger.warning("Output folder already exists! {path}", path=self.path)
        else:
            os.makedirs(self.path)

        t = time.time()
        for shot, dataset in self.results.items():
            cdf = os.path.join(self.path, f"{shot}.nc")
            logger.trace("Saving result: {cdf}", cdf=cdf)
            dataset.to_netcdf(cdf)
        logger.info(
            "Saved results in {took:.3f} s: {path}",
            took=time.time() - t,
            path=self.path,
        )
        return self.path


class SingleOutputSetting(DictOutputSetting):
    """
    Abstract class that outputs data as a single object/file.
    """

    def __init__(self, path: str | bool = True):
        """
        Initialize empty SingleOutputSetting.

        Parameters
        ----------
        path : str | bool, default = True
            The path for writing results to disk.
            If True, a unique temporary location will be used.
            If False, no results are written to disk.
        """
        super().__init__(path=False)
        self.result = None
        if path is True:
            ext = ".csv" if "DataFrame" in self.__class__.__name__ else ".nc"
            path = os.path.join(get_temporary_folder(), f"output{ext}")
            if os.path.exists(path) or "pytest" in sys.modules:
                _, path = tempfile.mkstemp(
                    dir=get_temporary_folder(), prefix="output-", suffix=ext
                )
        self.path = path

    @abstractmethod
    def concat(self) -> OutputSingleType:
        """
        Concatenate the resulting object.

        Returns
        -------
        xr.Dataset | xr.DataTree | pd.DataFrame
            The resulting object.
        """

    def get_results(self) -> OutputSingleType:
        """
        Get the resulting object.

        Returns
        -------
        xr.Dataset | xr.DataTree | pd.DataFrame
            The resulting object.
        """
        logger.debug("Concatenating {tot} shots.", tot=len(self.results))
        self.result = self.concat()
        self.results = {}
        return self.result

    def to_disk(self) -> str:
        """
        Save the resulting object into a file.
        """

        if not self.path:
            return ""
        if os.path.exists(self.path) and os.path.getsize(self.path):
            raise FileExistsError(f"File already exists! {self.path}")

        logger.debug(
            "Saving {type}: {path}", type=self.result.__class__.__name__, path=self.path
        )

        t = time.time()
        for method in ["to_netcdf", "to_csv"]:
            if not hasattr(self.result, method):
                continue
            getattr(self.result, method)(self.path)
            break
        else:
            raise NotImplementedError("Could not save object to file.")
        logger.info(
            "Saved {type} in {took:.3f} s: {path}",
            type=self.result.__class__.__name__,
            took=time.time() - t,
            path=self.path,
        )
        return self.path


class DatasetOutputSetting(SingleOutputSetting):
    """
    Outputs data as a single Dataset.
    """

    def concat(self) -> xr.Dataset:
        """
        Concatenate the resulting Dataset.

        Returns
        -------
        xr.Dataset
            The resulting Dataset.
        """
        return xr.concat(self.results.values(), dim="idx")


class DataTreeOutputSetting(SingleOutputSetting):
    """
    Outputs data as a single DataTree.
    """

    def concat(self) -> xr.DataTree:
        """
        Concatenate the resulting DataTree.

        Returns
        -------
        xr.DataTree
            The DataTree containing the results, with shots as keys.
        """
        return xr.DataTree.from_dict(self.results)


class DataFrameOutputSetting(SingleOutputSetting):
    """
    Outputs data as a DataFrame.
    """

    def concat(self) -> pd.DataFrame:
        """
        Concatenate the resulting DataFrame.

        Returns
        -------
        pd.DataFrame
            The resulting DataFrame.
        """
        base = ["shot", "time"]
        df = xr.concat(self.results.values(), dim="idx").to_dataframe()
        cols = base + [c for c in sorted(df.columns) if c not in base]
        return df[cols]


# --8<-- [start:output_setting_dict]
_output_setting_mappings: Dict[str, OutputSetting] = {
    "dataframe": DataFrameOutputSetting,
    "dataset": DatasetOutputSetting,
    "datatree": DataTreeOutputSetting,
    "df": DataFrameOutputSetting,
    "dict": DictOutputSetting,
    "ds": DatasetOutputSetting,
    "dt": DataTreeOutputSetting,
    "pandas": DataFrameOutputSetting,
    "pd": DataFrameOutputSetting,
    "xarray": DatasetOutputSetting,
    "xr": DatasetOutputSetting,
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
            return output_setting_object()
        # check suffixes
        for suffix, output_setting_type in _file_suffix_to_output_setting.items():
            if output_setting.endswith(suffix):
                return output_setting_type(output_setting)

    if isinstance(output_setting, list):
        return OutputSettingList(output_setting)

    raise ValueError(f"Invalid output processor {output_setting}")
