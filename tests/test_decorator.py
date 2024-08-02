import numpy as np
import pytest

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


@pytest.mark.parametrize("tokamak", [None, resolve_tokamak_from_environment()])
def test_tokamak_parameter(tokamak):
    """Ensure physics methods run when the tokamak parameter is set as either
    `None` or a specific tokamak."""
    col_name = "x"

    @physics_method(columns=[col_name], tokamak=tokamak)
    def my_physics_method(params: PhysicsMethodParams):
        return {col_name: np.ones(shape=len(params.times))}

    retrieval_settings = RetrievalSettings(
        run_tags=[],
        run_columns=[col_name],
        only_requested_columns=True,
        custom_physics_methods=[my_physics_method],
    )
    shot_data = get_shots_data(
        tokamak="cmod",
        shotlist_setting=[1150805012],
        retrieval_settings=retrieval_settings,
        output_setting="list",
        num_processes=1,
    )
    assert col_name in shot_data[0]
