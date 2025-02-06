#!/usr/bin/env python3

"""
Unit tests for ensuring data can be outputted in multiple formats including
lists, dictionaries, DataFrames, csv, hdf5, and to an SQL table.
"""

import os
from typing import Dict

import pandas as pd
import pytest

from disruption_py.config import config
from disruption_py.core.utils.misc import without_duplicates
from disruption_py.inout.sql import ShotDatabase
from disruption_py.settings.output_setting import (
    BatchedCSVOutputSetting,
    SQLOutputSetting,
)
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_database, get_shots_data
from tests.conftest import skip_on_fast_execution
from tests.utils.data_difference import assert_frame_equal_unordered

WRITE_DATABASE_TABLE_NAME = config().inout.sql.write_database_table_name
FIRST_ITERATION_COLUMNS = config().inout.sql.protected_columns + ["beta_p"]
SECOND_ITERATION_COLUMNS = ["kappa"]
ALL_ITERATION_COLUMNS = without_duplicates(
    FIRST_ITERATION_COLUMNS + SECOND_ITERATION_COLUMNS
)


@pytest.fixture(scope="module", name="shot_database")
def shot_database_fixture(tokamak) -> ShotDatabase:
    """
    Fixture for creating a ShotDatabase instance.
    """
    return get_database(tokamak=tokamak)


@pytest.fixture(autouse=True, scope="module")
def setup_shot_database(shotlist, shot_database):
    """
    This fixture automatically removes data for each shot in the
    provided shotlist from the shot_database.
    """
    for shot in shotlist:
        shot_database.remove_shot_data(shot)


@pytest.fixture(scope="module", name="initial_mdsplus_data")
def initial_mdsplus_data_fixture(shotlist, tokamak, test_file_path_f) -> Dict:
    """Get data in multiple formats"""
    output_settings = [
        "dataframe",
        "list",
        "dict",
        test_file_path_f(".csv"),
        test_file_path_f(".hdf5"),
        SQLOutputSetting(table_name=WRITE_DATABASE_TABLE_NAME),
    ]

    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        run_columns=FIRST_ITERATION_COLUMNS,
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
    df_output, list_output, dict_output, csv_output, hdf_output, sql_output = (
        initial_mdsplus_data
    )
    assert isinstance(df_output, pd.DataFrame), "DataFrame output does not exist"
    assert isinstance(list_output, list), "List output does not exist"
    assert isinstance(dict_output, dict), "Dict output does not exist"
    assert isinstance(csv_output, pd.DataFrame), "DataFrame from CSV does not exist"
    assert isinstance(hdf_output, pd.DataFrame), "DataFrame from HDF5 does not exist"
    assert isinstance(sql_output, pd.DataFrame), "DataFrame from SQL does not exist"
    assert os.path.exists(test_file_path_f(".csv")), ".csv output does not exist"
    assert os.path.exists(test_file_path_f(".hdf5")), ".hdf5 output does not exist"

    assert_frame_equal_unordered(df_output, csv_output)
    assert_frame_equal_unordered(csv_output, hdf_output)
    assert_frame_equal_unordered(hdf_output, sql_output)


@skip_on_fast_execution
def test_sql_output_setting(
    shotlist, shot_database: ShotDatabase, initial_mdsplus_data_df, test_file_path_f
) -> Dict:
    """
    Ensure SQL output setting works by reading from an initial writeback to
    SQL, then by updating the SQL with another writeback, and making sure data from
    MDSplus matches the data from SQL for both retrievals.
    """

    # Test initial database writeback
    result = shot_database.get_shots_data(
        shotlist=shotlist,
        cols=FIRST_ITERATION_COLUMNS,
        sql_table=WRITE_DATABASE_TABLE_NAME,
    )
    try:
        assert_frame_equal_unordered(
            result[FIRST_ITERATION_COLUMNS],
            initial_mdsplus_data_df[FIRST_ITERATION_COLUMNS],
        )
    except AssertionError:
        print("Initial writeback to SQL failed!")
        result.to_csv(test_file_path_f("-1L.csv"))
        initial_mdsplus_data_df[FIRST_ITERATION_COLUMNS].to_csv(
            test_file_path_f("-1R.csv")
        )
        raise

    # Do second retrieval that updates the data for the columns
    retrieval_settings = RetrievalSettings(
        efit_nickname_setting="disruption",
        run_columns=ALL_ITERATION_COLUMNS,
        only_requested_columns=True,
    )
    shot_data = get_shots_data(
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=SQLOutputSetting(
            table_name=WRITE_DATABASE_TABLE_NAME,
            should_override_columns=SECOND_ITERATION_COLUMNS,
        ),
        num_processes=2,
    )
    result = shot_database.get_shots_data(
        shotlist=shotlist,
        cols=ALL_ITERATION_COLUMNS,
        sql_table=WRITE_DATABASE_TABLE_NAME,
    )
    try:
        assert_frame_equal_unordered(
            result[ALL_ITERATION_COLUMNS], shot_data[ALL_ITERATION_COLUMNS]
        )
    except AssertionError:
        print("Second writeback to SQL failed!")
        result[ALL_ITERATION_COLUMNS].to_csv(test_file_path_f("-2L.csv"))
        shot_data[ALL_ITERATION_COLUMNS].to_csv(test_file_path_f("-2R.csv"))
        raise


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
