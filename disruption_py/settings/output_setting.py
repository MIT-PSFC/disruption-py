#!/usr/bin/env python3

"""
Handles output settings for retrieving and saving shot data.

This module provides classes and methods to manage various output settings.
"""

import os
import tempfile
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Type, TypeAlias, Union

import pandas as pd
import xarray as xr
from loguru import logger

from disruption_py.core.utils.misc import get_max_rss, get_temporary_folder, shot_msg
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


OutputSettingType: TypeAlias = Union["OutputSetting", str, List["OutputSettingType"]]
OutputDictType: TypeAlias = Dict[int, xr.Dataset]
OutputSingleType: TypeAlias = xr.Dataset | xr.DataTree | pd.DataFrame
OutputType: TypeAlias = OutputDictType | OutputSingleType


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
    def to_disk(self) -> str | List[str]:
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

        # include DictOutputSetting but exclude subclasses
        # pylint: disable=unidiomatic-typecheck

        self.results: Dict[int, xr.Dataset] = {}
        self.shards: Dict[int, str] = {}

        if path is True:
            path = os.path.join(get_temporary_folder(), "output")

        if type(self) is not DictOutputSetting:
            # subclasses write/remove shards in the temporary folder
            pass
        elif path is False:
            # explicit DictOutputSetting without path defeats its purpose
            logger.warning("Memory optimization requires unloading to disk!")
        else:
            # standard
            logger.trace("Creating output folder: {path}", path=path)
            os.makedirs(path, exist_ok=True)
            if os.listdir(path):
                logger.warning("Output folder is not empty! {path}", path=path)

        self.path = path

    def _output_shot(self, params: OutputSettingParams):
        """
        Store a single result in the dictionary.

        Parameters
        ----------
        params : OutputSettingParams
            The parameters for outputting shot results.
        """

        # include DictOutputSetting but exclude subclasses
        # pylint: disable=unidiomatic-typecheck

        if type(self) is DictOutputSetting and self.path is False:
            # do not shard
            self.results[params.shot_id] = params.result
            return

        file = f"{params.shot_id}.nc"
        if type(self) is DictOutputSetting:
            # shard (i.e. store) into output folder
            shard = os.path.join(self.path, file)
            if os.path.exists(shard):
                logger.warning(f"Output file already exists! {shard}")
                # rename shard to avoid losing data
                fd, shard = tempfile.mkstemp(
                    dir=self.path, prefix=f"{params.shot_id}.", suffix=".nc"
                )
                os.close(fd)
        else:
            # shard into temporary folder
            fd, shard = tempfile.mkstemp(
                dir=get_temporary_folder(), prefix=f".{params.shot_id}.", suffix=".nc"
            )
            os.close(fd)

        # save to disk
        self.shards[params.shot_id] = shard
        logger.trace(
            shot_msg("Saving shard: {shard}"), shot=params.shot_id, shard=shard
        )
        params.result.to_netcdf(shard)

        # lazy reload
        params.result.close()
        self.results[params.shot_id] = params.result = xr.open_dataset(shard)

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
        # if the DictOutputSetting.to_disk method is called,
        # there is nothing to do but to log the output folder
        # subclasses will handle their own writing to disk
        logger.info("Saved results: {path}", path=self.path)
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
            If False, no results are written to disk (excluding temporary shards).
        """

        # although we instantiate a DictOutputSetting without path, temporary
        # shards are still saved to disk to protect against data loss, and cleaned
        # up only after successfully persisting the SingleOutputSetting to disk
        super().__init__(path=False)
        self.result = None

        if path is True:
            ext = "csv" if isinstance(self, DataFrameOutputSetting) else "nc"
            path = os.path.join(get_temporary_folder(), f"output.{ext}")

        if path and os.path.exists(path):
            logger.warning(f"Output file already exists! {path}")
            # rename file to avoid losing data
            folder = os.path.dirname(path)
            name, ext = os.path.splitext(os.path.basename(path))
            _, path = tempfile.mkstemp(dir=folder, prefix=f"{name}.", suffix=ext)

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

        logger.debug("Reading {tot:,} shots...", tot=len(self.results))
        took = -time.time()
        for result in self.results.values():
            result.load()
        took += time.time()
        logger.info(
            "Read {tot:,} shots in {sec:.3f}s.", tot=len(self.results), sec=took
        )
        logger.debug(
            "Read shots: MaxRSS = {mem:,.1f} MB",
            mem=get_max_rss(),
        )

        logger.debug("Concatenating {tot:,} shots...", tot=len(self.results))
        took = -time.time()
        self.result = self.concat()
        took += time.time()
        logger.info(
            "Concatenated {tot:,} shots in {sec:.3f}s.", tot=len(self.results), sec=took
        )
        logger.debug(
            "Concatenated shots: MaxRSS = {mem:,.1f} MB",
            mem=get_max_rss(),
        )
        self.results = {}

        return self.result

    def to_disk(self) -> str:
        """
        Save the resulting object into a file.
        """

        if self.path:
            os.makedirs(os.path.dirname(self.path), exist_ok=True)
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
                type=type(self).__name__.removesuffix("OutputSetting"),
                took=time.time() - t,
                path=self.path,
            )

        for shard in self.shards.values():
            logger.trace("Removing shard: {shard}", shard=shard)
            os.remove(shard)

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
        if not self.results:
            logger.critical("Nothing to concatenate!")
            return xr.Dataset()

        ds = xr.concat(self.results.values(), dim="idx", combine_attrs="no_conflicts")
        if "shot" not in ds.coords or "time" not in ds.coords:
            return ds

        took = -time.time()
        ds = ds.sortby(["shot", "time"])
        took += time.time()
        logger.debug("Sorted {tot:,} rows in {sec:.3f}s.", tot=len(ds.idx), sec=took)
        return ds


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
        if not self.results:
            logger.critical("Nothing to concatenate!")
            return xr.DataTree()
        return xr.DataTree.from_dict({str(k): v for k, v in self.results.items()})


class DataFrameOutputSetting(DatasetOutputSetting):
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
        if not self.results:
            logger.critical("Nothing to concatenate!")
            return pd.DataFrame()
        df = super().concat().to_dataframe()
        base = ["shot", "time"]
        cols = base + [c for c in sorted(df.columns) if c not in base]
        return df[cols]


# --8<-- [start:output_setting_dict]
_output_setting_mappings: Dict[str, Type[OutputSetting]] = {
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
