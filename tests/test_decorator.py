#!/usr/bin/env python3

"""
This module contains tests to ensure that physics methods can be executed
correctly when the tokamak parameter is set to either `None` or a specific
tokamak instance.
"""

import os

import pytest

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data


# Module-scope so spawn/forkserver workers can unpickle the function reference.
@physics_method(columns=["x"])
def _physics_method_any_tokamak(params: PhysicsMethodParams):
    return {"x": params.times**0}


@physics_method(columns=["x"], tokamak=resolve_tokamak_from_environment())
def _physics_method_specific_tokamak(params: PhysicsMethodParams):
    return {"x": params.times**0}


@pytest.mark.parametrize(
    "method",
    [_physics_method_any_tokamak, _physics_method_specific_tokamak],
    ids=["tokamak=None", "tokamak=resolved"],
)
def test_tokamak_parameter(shotlist, method, test_folder_f):
    """
    Ensure physics methods run when the tokamak parameter is set as either
    `None` or a specific tokamak.
    """
    retrieval_settings = RetrievalSettings(
        run_columns=["x"],
        only_requested_columns=True,
        custom_physics_methods=[method],
    )
    shot_data = get_shots_data(
        shotlist_setting=shotlist[:1],
        retrieval_settings=retrieval_settings,
        output_setting=os.path.join(test_folder_f, "output.nc"),
        log_settings=LogSettings(
            console_level="WARNING",
            file_path=os.path.join(test_folder_f, "output.log"),
        ),
    )
    assert "x" in shot_data.data_vars
    assert shotlist[0] in shot_data.shot
    assert len(shot_data.time)
    assert all(shot_data["x"].values == 1)
