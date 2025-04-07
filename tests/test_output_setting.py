#!/usr/bin/env python3

"""
Unit tests for ensuring data can be outputted in multiple formats including
lists, dictionaries, DataFrames, csv, hdf5, and to an SQL table.
"""

import os
from typing import Dict

import pandas as pd
import pytest

from disruption_py.settings.output_setting import BatchedCSVOutputSetting
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.data_difference import assert_frame_equal_unordered


@pytest.fixture(scope="module", name="initial_mdsplus_data")
def initial_mdsplus_data_fixture(shotlist, tokamak, test_file_path_f) -> Dict:
    """Get data in multiple formats"""
    output_settings = [
        "dataframe",
        "dict",
        test_file_path_f(".csv"),
    ]

    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        run_columns=["kappa"],
        only_requested_columns=True,
    )
    all_outputs = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=output_settings,
        num_processes=2,
    )
    return all_outputs


@pytest.fixture(scope="module", name="initial_mdsplus_data_df")
def initial_mdsplus_data_df_fixture(initial_mdsplus_data):
    """Return only the dataframe output type"""
    return initial_mdsplus_data[0]


def test_output_exists(initial_mdsplus_data, test_file_path_f):
    """
    Test creation of all output formats except SQL.
    """
    df_output, dict_output, csv_output = initial_mdsplus_data
    assert isinstance(df_output, pd.DataFrame), "DataFrame output does not exist"
    assert isinstance(dict_output, dict), "Dict output does not exist"
    assert isinstance(csv_output, pd.DataFrame), "DataFrame from CSV does not exist"
    assert os.path.exists(test_file_path_f(".csv")), ".csv output does not exist"

    assert_frame_equal_unordered(df_output, csv_output)


def test_batch_csv(tokamak, test_file_path_f, shotlist):
    """
    Test the batch csv output setting to ensure it outputs the same columns in
    the same order as the dataframe in memory.
    """
    csv = test_file_path_f("-batch.csv")
    out = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        num_processes=2,
        # Use a batch size less than the number of shots to ensure multiple batches
        # are written to the CSV file.
        output_setting=BatchedCSVOutputSetting(filepath=csv, batch_size=1),
    )
    df = pd.read_csv(csv)
    assert_frame_equal_unordered(out, df)
