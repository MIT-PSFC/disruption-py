import os
from typing import Dict

import pandas as pd
import pytest

from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.output_type_request import SQLOutputRequest
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.utils.constants import BASE_PROTECTED_COLUMNS
from disruption_py.utils.utils import without_duplicates

TEST_SHOTS = [
    1150805012,  # Flattop Disruption
    1150805013,  # No Disruption
]

TABLE_NAME = "disruption_warning_test"
FIRST_ITERATION_COLUMNS = BASE_PROTECTED_COLUMNS + ["beta_p"]
SECOND_ITERATION_COLUMNS = ["kappa"]
ALL_ITERATION_COLUMNS = without_duplicates(
    FIRST_ITERATION_COLUMNS + SECOND_ITERATION_COLUMNS
)

@pytest.fixture(scope="class")
def cmod_handler():
    return CModHandler()

@pytest.fixture(scope="class")
def initial_mdsplus_data(cmod_handler: CModHandler, shotlist) -> Dict:
    shot_settings = ShotSettings(
        set_times_request="efit",
        efit_tree_name="efit18",
        run_columns=FIRST_ITERATION_COLUMNS,
        only_requested_columns=True,
    )
    shot_data, _ = cmod_handler.get_shots_data(
        shot_ids_request=shotlist,
        shot_settings=shot_settings,
        output_type_request=["dataframe", SQLOutputRequest(table_name=TABLE_NAME)],
        num_processes=1,
    )
    return shot_data


def assert_frame_equal_unordered(df1: pd.DataFrame, df2: pd.DataFrame):
    df1_sorted = df1.sort_values(by=BASE_PROTECTED_COLUMNS).reset_index(drop=True)
    df2_sorted = df2.sort_values(by=BASE_PROTECTED_COLUMNS).reset_index(drop=True)
    print(df1_sorted)
    print(df2_sorted)
    pd.testing.assert_frame_equal(df1_sorted, df2_sorted, check_like=True)


@pytest.mark.skipif(
    os.path.exists("/fusion/projects/disruption_warning"), reason="on DIII-D"
)
def test_update_data(cmod_handler: CModHandler, shotlist, initial_mdsplus_data) -> Dict:
    # Test initial database readback
    cmod_database = cmod_handler.database
    result = cmod_database.get_shots_data(shotlist, sql_table=TABLE_NAME)
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
    shot_data, _ = cmod_handler.get_shots_data(
        shot_ids_request=shotlist,
        shot_settings=shot_settings,
        output_type_request=[
            "dataframe",
            SQLOutputRequest(
                table_name=TABLE_NAME, should_override_columns=SECOND_ITERATION_COLUMNS
            ),
        ],
        num_processes=1,
    )
    cmod_database = cmod_handler.database
    result = cmod_database.get_shots_data(shotlist, sql_table=TABLE_NAME)
    assert_frame_equal_unordered(
        result[ALL_ITERATION_COLUMNS], shot_data[ALL_ITERATION_COLUMNS]
    )
