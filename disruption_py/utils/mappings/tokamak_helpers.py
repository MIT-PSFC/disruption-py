#!/usr/bin/env python3

from logging import Logger
import os
from typing import Callable

from disruption_py.shots.cmod_shot_manager import CModShotManager
from disruption_py.shots.d3d_shot_manager import D3DShotManager
from disruption_py.utils.constants import (
    EXPECTED_FAILURE_COLUMNS,
    TEST_COLUMNS,
    TEST_SHOTS,
)
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
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


def get_tokamak_from_environment():
    if "DISPY_TOKAMAK" in os.environ:
        return Tokamak[os.environ["DISPY_TOKAMAK"]]
    if os.path.exists("/usr/local/mfe/disruptions"):
        return Tokamak.CMOD
    if os.path.exists("/fusion/projects/disruption_warning"):
        return Tokamak.D3D
    return None


def resolve_tokamak(tokamak: Tokamak, logger: Logger = None):
    if tokamak is None:
        tokamak = get_tokamak_from_environment()
        if logger:
            logger.info(f"No tokamak argument given. Detected tokamak: {tokamak.value}")
        return tokamak
    else:
        return map_string_to_enum(tokamak, Tokamak)


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
