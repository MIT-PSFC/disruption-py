import os
from typing import Dict

import pandas as pd
import pytest

from disruption_py.database import ShotDatabase
from disruption_py.main import get_database, get_shots_data
from disruption_py.settings.output_type_request import SQLOutputRequest
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.utils.constants import (
    BASE_PROTECTED_COLUMNS,
)
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.utils import without_duplicates

WRITE_DATABASE_TABLE_NAME = "disruption_warning_test"  # overwrite value in constants
FIRST_ITERATION_COLUMNS = BASE_PROTECTED_COLUMNS + ["beta_p"]
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
    shot_settings = ShotSettings(
        set_times_request="efit",
        efit_tree_name="efit18",
        run_columns=FIRST_ITERATION_COLUMNS,
        only_requested_columns=True,
    )
    shot_data, _ = get_shots_data(
        tokamak=tokamak,
        shot_ids_request=shotlist,
        shot_settings=shot_settings,
        output_type_request=[
            "dataframe",
            SQLOutputRequest(table_name=WRITE_DATABASE_TABLE_NAME),
        ],
        num_processes=1,
    )
    return shot_data


def assert_frame_equal_unordered(df1: pd.DataFrame, df2: pd.DataFrame):
    df1_sorted = df1.sort_values(by=BASE_PROTECTED_COLUMNS).reset_index(drop=True)
    df2_sorted = df2.sort_values(by=BASE_PROTECTED_COLUMNS).reset_index(drop=True)
    print(df1_sorted)
    print(df2_sorted)
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

    # do second request that updates the data for the columns
    shot_settings = ShotSettings(
        set_times_request="efit",
        efit_tree_name="efit18",
        run_columns=ALL_ITERATION_COLUMNS,
        only_requested_columns=True,
    )
    shot_data, _ = get_shots_data(
        shot_ids_request=shotlist,
        shot_settings=shot_settings,
        output_type_request=[
            "dataframe",
            SQLOutputRequest(
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
