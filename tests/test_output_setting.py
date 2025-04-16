#!/usr/bin/env python3

"""
Unit tests for ensuring data can be outputted in multiple formats including
lists, dictionaries, DataFrames, csv, hdf5, and to an SQL table.
"""

import os
from typing import Dict

import pandas as pd
import pytest
import xarray as xr

from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.output_setting import DataTreeOutputSetting
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


@pytest.fixture(scope="module", name="fresh_data")
def fresh_data_fixture(shotlist, tokamak, test_folder_m) -> Dict:
    """
    Get data in multiple formats.
    """
    output_settings = [
        os.path.join(test_folder_m, "output/"),
        os.path.join(test_folder_m, "dataset.nc"),
        os.path.join(test_folder_m, "dataframe.csv"),
        DataTreeOutputSetting(path=os.path.join(test_folder_m, "datatree.nc")),
    ]
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        run_columns=["kappa"],
        only_requested_columns=True,
    )
    return get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=output_settings,
        log_settings=LogSettings(
            console_level="WARNING",
            log_file_path=os.path.join(test_folder_m, "output.log"),
        ),
        num_processes=2,
    )


def test_output_exists(fresh_data, test_folder_m):
    """
    Test creation of all output formats.
    """

    # run and get from memory
    dict_out, ds_out, df_out, dt_out = fresh_data

    # get from disk
    dict_path = os.path.join(test_folder_m, "output/")
    ds_path = os.path.join(test_folder_m, "dataset.nc")
    dt_path = os.path.join(test_folder_m, "datatree.nc")
    df_path = os.path.join(test_folder_m, "dataframe.csv")
    ds_dsk = xr.open_dataset(ds_path)
    dt_dsk = xr.open_datatree(dt_path)
    df_dsk = pd.read_csv(df_path, index_col=0)

    # path existence
    assert os.path.exists(dict_path), "Could not find dict folder"
    assert os.path.exists(ds_path), "Could not find dataset file"
    assert os.path.exists(dt_path), "Could not find datatree file"
    assert os.path.exists(df_path), "Could not find dataframe file"

    # format types
    assert isinstance(dict_out, dict), "Wrong type for dict output"
    assert isinstance(ds_out, xr.Dataset), "Wrong type for Dataset output"
    assert isinstance(ds_dsk, xr.Dataset), "Wrong type for Dataset output"
    assert isinstance(dt_out, xr.DataTree), "Wrong type for DataTree output"
    assert isinstance(dt_dsk, xr.DataTree), "Wrong type for DataTree output"
    assert isinstance(df_out, pd.DataFrame), "Wrong type for DataFrame output"
    assert isinstance(df_dsk, pd.DataFrame), "Wrong type for DataFrame output"

    # disk equivalence
    xr.testing.assert_identical(ds_out, ds_dsk)
    xr.testing.assert_identical(dt_out, dt_dsk)
    pd.testing.assert_frame_equal(df_out, df_dsk)

    # format equivalence
    xr.testing.assert_identical(ds_out, xr.concat(dict_out.values(), dim="idx"))
    xr.testing.assert_identical(
        ds_out, xr.concat([dt.to_dataset() for dt in dt_out.values()], dim="idx")
    )
    pd.testing.assert_frame_equal(df_out, ds_out.to_dataframe()[df_out.columns])
