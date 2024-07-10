#!/usr/bin/env python3

import os
from typing import Dict

import pandas as pd
import pytest

from disruption_py.config import config
from disruption_py.io.sql import ShotDatabase
from disruption_py.workflow import get_database, get_shots_data
from disruption_py.settings.output_setting import SQLOutputSetting
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.machine.tokamak import Tokamak
from disruption_py.core.utils.misc import without_duplicates

WRITE_DATABASE_TABLE_NAME = config().database.write_database_table_name
FIRST_ITERATION_COLUMNS = config().database.protected_columns + ["beta_p"]
SECOND_ITERATION_COLUMNS = ["kappa"]
ALL_ITERATION_COLUMNS = without_duplicates(
    FIRST_ITERATION_COLUMNS + SECOND_ITERATION_COLUMNS
)


@pytest.fixture(scope="class")
def shot_database(tokamak) -> ShotDatabase:
    return get_database(tokamak=tokamak)


@pytest.fixture(autouse=True, scope="class")
def setup_shot_database(shotlist, shot_database):
    for shot in shotlist:
        shot_database.remove_shot_data(shot)


@pytest.fixture(scope="class")
def initial_mdsplus_data(shotlist, tokamak, setup_shot_database) -> Dict:
    if tokamak is Tokamak.D3D:
        pytest.skip("Skipping test on DIII-D")
    retrieval_settings = RetrievalSettings(
        time_setting="efit",
        efit_tree_name="efit18",
        run_columns=FIRST_ITERATION_COLUMNS,
        only_requested_columns=True,
    )
    shot_data, _ = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=[
            "dataframe",
            SQLOutputSetting(table_name=WRITE_DATABASE_TABLE_NAME),
        ],
        num_processes=1,
    )
    return shot_data


def assert_frame_equal_unordered(df1: pd.DataFrame, df2: pd.DataFrame):
    df1_sorted = df1.sort_values(by=config().database.protected_columns).reset_index(
        drop=True
    )
    df2_sorted = df2.sort_values(by=config().database.protected_columns).reset_index(
        drop=True
    )
    pd.testing.assert_frame_equal(df1_sorted, df2_sorted, check_like=True)


def test_update_data(
    shotlist, initial_mdsplus_data, tokamak, shot_database: ShotDatabase
) -> Dict:
    if tokamak is Tokamak.D3D:
        pytest.skip("Skipping test on DIII-D")
    # Test initial database readback
    result = shot_database.get_shots_data(shotlist, sql_table=WRITE_DATABASE_TABLE_NAME)
    assert_frame_equal_unordered(
        result[FIRST_ITERATION_COLUMNS], initial_mdsplus_data[FIRST_ITERATION_COLUMNS]
    )

    # do second retrieval that updates the data for the columns
    retrieval_settings = RetrievalSettings(
        time_setting="efit",
        efit_tree_name="efit18",
        run_columns=ALL_ITERATION_COLUMNS,
        only_requested_columns=True,
    )
    shot_data, _ = get_shots_data(
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        output_setting=[
            "dataframe",
            SQLOutputSetting(
                table_name=WRITE_DATABASE_TABLE_NAME,
                should_override_columns=SECOND_ITERATION_COLUMNS,
            ),
        ],
        num_processes=1,
    )
    result = shot_database.get_shots_data(shotlist, sql_table=WRITE_DATABASE_TABLE_NAME)
    assert_frame_equal_unordered(
        result[ALL_ITERATION_COLUMNS], shot_data[ALL_ITERATION_COLUMNS]
    )
