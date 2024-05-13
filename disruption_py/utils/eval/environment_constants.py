import os
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.handlers.d3d_handler import D3DHandler
from disruption_py.utils.constants import CMOD_EXPECTED_FAILURE_COLUMNS, CMOD_TEST_COLUMNS, CMOD_TEST_SHOTS, D3D_EXPECTED_FAILURE_COLUMNS, D3D_TEST_COLUMNS, D3D_TEST_SHOTS
from disruption_py.utils.mappings.tokamak import Tokamak


def get_test_expected_failure_columns(tokamak : Tokamak):
    if tokamak == Tokamak.CMOD:
        return CMOD_EXPECTED_FAILURE_COLUMNS
    elif tokamak == Tokamak.D3D:
        return D3D_EXPECTED_FAILURE_COLUMNS
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))


def get_test_handler(tokamak : Tokamak):
    if tokamak is Tokamak.CMOD:
        return CModHandler()
    elif tokamak is Tokamak.D3D:
        return D3DHandler()
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))


def get_test_shot_ids(tokamak : Tokamak) -> list[int]:
    if tokamak == Tokamak.CMOD:
        shot_id_dict = CMOD_TEST_SHOTS
    elif tokamak == Tokamak.D3D:
        shot_id_dict =  D3D_TEST_SHOTS
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))
    
    if "GITHUB_ACTIONS" in os.environ:
        shot_id_dict = {key: value for key, value in shot_id_dict.items() if "_fast" in key}
        
    return list(shot_id_dict.values())


def get_test_columns(tokamak : Tokamak):
    if tokamak == Tokamak.CMOD:
        return CMOD_TEST_COLUMNS
    elif tokamak == Tokamak.D3D:
        return D3D_TEST_COLUMNS
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))