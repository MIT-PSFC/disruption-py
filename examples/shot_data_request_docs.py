#!/usr/bin/env python3

""" Used in the documentation for the physics methods. """

import numpy as np
import pandas as pd

from disruption_py.core.physics_method.decorator import physics_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak


@physics_method(columns=["upper_gap", "lower_gap"], tokamak=Tokamak.D3D)
def decorated_physics_method(self, params: PhysicsMethodParams) -> pd.DataFrame:
    """All parametered methods passed to `get_shots_data` will be called once for every shot retrieved.
    Decorated methods may call other decorated methods, however, execution order is not guranteed as calls
    will be reordered to minimize resource usage based on the `physics_method` decorator.

    Parameters
    ----------
    params : PhysicsMethodParams
        Parameters passed by disruption_py to the decorated method that should be used to help retrieve the shot data from MDSplus.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the results of the decorated method, with each returned parameter being a column.
        The dataframe should contain the same number of rows as the timebase (`physics_method_params.shot.times`).
    """
    pass


# Paramater cached method example
# --8<-- [start:kappa_area_request_example]


class KappaAreaRequest:

    @staticmethod
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

        aminor[aminor <= 0] = 0.001  # make sure aminor is not 0 or less than 0
        # make sure area is not 0 or less than 0
        area[area <= 0] = 3.14 * 0.001**2
        return pd.DataFrame(
            {
                "kappa_area": params.interpolation_method(
                    times, area / (np.pi * aminor**2), params.shot_props.times
                )
            }
        )


# --8<-- [end:kappa_area_request_example]
