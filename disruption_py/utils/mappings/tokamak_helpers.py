#!/usr/bin/env python3

import os

from disruption_py.shots.cmod_shot_manager import CModShotManager
from disruption_py.shots.d3d_shot_manager import D3DShotManager
from disruption_py.utils.constants import (
    EXPECTED_FAILURE_COLUMNS,
    TEST_COLUMNS,
    TEST_SHOTS,
)
from disruption_py.utils.mappings.tokamak import Tokamak


def get_tokamak_from_shot_id(shot_id):
    if isinstance(shot_id, str):
        shot_len = len(shot_id)
    elif isinstance(shot_id, int):
        # math.log10 is faster and safer for large numbers but we assume shot_id is relatively small
        shot_len = len(str(shot_id))
    else:
        raise ValueError(f"shot_id must be a string or integer, not {type(shot_id)}")

    if shot_len == 6:
        return Tokamak.D3D
    elif shot_len == 10:
        return Tokamak.CMOD
    else:
        raise NotImplementedError(f"Unable to handle shot_id of length {shot_len}")


def get_tokamak_shot_manager(tokamak: Tokamak):
    if tokamak == Tokamak.CMOD:
        return CModShotManager
    elif tokamak == Tokamak.D3D:
        return D3DShotManager
    else:
        raise ValueError("No shot manager for tokamak {}".format(tokamak))


def get_tokamak_test_expected_failure_columns(tokamak: Tokamak):
    return EXPECTED_FAILURE_COLUMNS.get(tokamak.value)


def get_tokamak_test_shot_ids(tokamak: Tokamak) -> list[int]:
    shot_id_dict = TEST_SHOTS.get(tokamak.value)

    if "GITHUB_ACTIONS" in os.environ:
        shot_id_dict = {
            key: value for key, value in shot_id_dict.items() if "_fast" in key
        }

    return list(shot_id_dict.values())


def get_tokamak_test_columns(tokamak: Tokamak):
    return TEST_COLUMNS.get(tokamak.value)


def built_in_method_factory(tokamak: Tokamak):
    if tokamak is Tokamak.D3D:
        from disruption_py.shots.parameter_methods.d3d.built_in import (
            D3D_DEFAULT_SHOT_DATA_REQUESTS,
        )

        return D3D_DEFAULT_SHOT_DATA_REQUESTS
    elif tokamak is Tokamak.CMOD:
        from disruption_py.shots.parameter_methods.cmod.built_in import (
            CMOD_DEFAULT_SHOT_DATA_REQUESTS,
        )

        return CMOD_DEFAULT_SHOT_DATA_REQUESTS
    else:
        raise ValueError(f"Invalid tokamak for built-ins {tokamak}")
