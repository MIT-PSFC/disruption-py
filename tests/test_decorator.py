#!/usr/bin/env python3

"""
This module contains tests to ensure that physics methods can be executed
correctly when the tokamak parameter is set to either `None` or a specific
tokamak instance.
"""

import os

import numpy as np
import pytest

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data


@pytest.mark.parametrize("tok", [None, resolve_tokamak_from_environment()])
def test_tokamak_parameter(shotlist, tok, test_folder_f):
    """
    Ensure physics methods run when the tokamak parameter is set as either
    `None` or a specific tokamak.
    """
    col_name = "x"

    # The physics method needs to be defined in the global scope because
    # multiprocessing and pickling don't work with locally-defined functions.
    # pylint: disable-next=global-variable-undefined
    global my_physics_method

    @physics_method(columns=[col_name], tokamak=tok)
    def my_physics_method(params: PhysicsMethodParams):
        return {col_name: np.ones(shape=len(params.times))}

    retrieval_settings = RetrievalSettings(
        run_columns=[col_name],
        only_requested_columns=True,
        custom_physics_methods=[my_physics_method],
    )
    shot_data = get_shots_data(
        shotlist_setting=shotlist[:1],
        retrieval_settings=retrieval_settings,
        output_setting=os.path.join(test_folder_f, "output.nc"),
        log_settings=LogSettings(
            console_log_level="WARNING",
            log_file_path=os.path.join(test_folder_f, "output.log"),
        ),
    )
    assert col_name in shot_data.data_vars
