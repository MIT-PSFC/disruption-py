#!/usr/bin/env python3

""" Used in the documentation for the physics methods. """

import numpy as np

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


@physics_method(columns=["upper_gap", "lower_gap"], tokamak=Tokamak.D3D)
def decorated_physics_method(params: PhysicsMethodParams) -> dict:
    """All parametrized methods passed to `get_shots_data` will be called once for every shot retrieved.
    Decorated methods may call other decorated methods, however, execution order is not guranteed as calls
    will be reordered to minimize resource usage based on the `physics_method` decorator.

    Parameters
    ----------
    params : PhysicsMethodParams
        Parameters passed by disruption_py to the decorated method that should be used to help retrieve the shot data from MDSplus.

    Returns
    -------
    dict
        Dictionary containing the results of the decorated method, with each returned
        parameter being a key-value pair. Each of the dictionary's values should be the same
        length as the timebase (`params.times`).
    """
    pass


# # Paramater cached method example
# # --8<-- [start:kappa_area_request_example]


# pylint: disable=duplicate-code
@physics_method(columns=["kappa_area"], tokamak=Tokamak.CMOD)
def _get_kappa_area(params: PhysicsMethodParams):
    aminor = params.mds_conn.get_data(
        r"\efit_aeqdsk:aminor", tree_name="_efit_tree", astype="float64"
    )
    area = params.mds_conn.get_data(
        r"\efit_aeqdsk:area", tree_name="_efit_tree", astype="float64"
    )
    times = params.mds_conn.get_data(
        r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
    )

    # Ensure aminor and area are not 0 or less than 0
    aminor[aminor <= 0] = 0.001
    area[area <= 0] = 3.14 * 0.001**2
    return {
        "kappa_area": params.interpolation_method(
            times, area / (np.pi * aminor**2), params.times
        )
    }


# --8<-- [end:kappa_area_request_example]


retrieval_settings = RetrievalSettings(
    time_setting="efit",
    # run only the custom kappa_area method
    run_tags=[],
    run_columns=["kappa_area"],
    only_requested_columns=True,
    custom_physics_methods=[_get_kappa_area],
)
shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=[1150805012],
    retrieval_settings=retrieval_settings,
    output_setting="list",
    num_processes=1,
)
