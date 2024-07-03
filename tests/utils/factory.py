import os
from disruption_py.machine.tokamak import Tokamak
from disruption_py.utils.constants import (
    EXPECTED_FAILURE_COLUMNS,
    TEST_COLUMNS,
    TEST_SHOTS,
)


def get_tokamak_test_expected_failure_columns(tokamak: Tokamak):
    return EXPECTED_FAILURE_COLUMNS.get(tokamak.value)


def get_tokamak_test_shotlist(tokamak: Tokamak) -> list[int]:
    shot_id_dict = TEST_SHOTS.get(tokamak.value)

    if "GITHUB_ACTIONS" in os.environ:
        shot_id_dict = {
            key: value for key, value in shot_id_dict.items() if "_fast" in key
        }

    return list(shot_id_dict.values())


def get_tokamak_test_columns(tokamak: Tokamak):
    return TEST_COLUMNS.get(tokamak.value)
