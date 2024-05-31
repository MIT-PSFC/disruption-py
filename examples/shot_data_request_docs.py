#!/usr/bin/env python3

""" Used in the documentation for the shot data request. """

import numpy as np
import pandas as pd

from disruption_py.settings.shot_data_request import (
    ShotDataRequest,
    ShotDataRequestParams,
)
from disruption_py.shots.helpers.method_caching import parameter_cached_method
from disruption_py.utils.mappings.tokamak import Tokamak


@parameter_cached_method(used_trees=["tree_1", "tree_2"])
def decorated_shot_data_method(self, params: ShotDataRequestParams) -> pd.DataFrame:
    """All decorated methods in subclasses of `ShotDataRequest` will be called once for every shot retrieved.
    Decorated methods may call other decorated methods, however, execution order is not guranteed as calls
    will be reordered to minimize resource usage based on the `parameter_cached_method` decorator.

    Parameters
    ----------
    params : ShotDataRequest
        Parameters passed by disruption_py to the decorated method that should be used to help retrieve the shot data from MDSplus.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the results of the decorated method, with each returned parameter being a column.
        The dataframe should contain the same number of rows as the timebase (`shot_data_request_params.shot.times`).
    """
    pass


# Paramater cached method example
# --8<-- [start:kappa_area_request_example]


class KappaAreaRequest(ShotDataRequest):

    @staticmethod
    @parameter_cached_method(
        columns=["kappa_area"], used_trees=["_efit_tree"], tokamak=Tokamak.CMOD
    )
    def _get_kappa_area(params: ShotDataRequestParams):
        aminor = params.shot_props.mds_conn.get_data(
            r"\efit_aeqdsk:aminor", tree_name="_efit_tree", astype="float64"
        )
        area = params.shot_props.mds_conn.get_data(
            r"\efit_aeqdsk:area", tree_name="_efit_tree", astype="float64"
        )
        times = params.shot_props.mds_conn.get_data(
            r"\efit_aeqdsk:time", tree_name="_efit_tree", astype="float64"
        )

        aminor[aminor <= 0] = 0.001  # make sure aminor is not 0 or less than 0
        # make sure area is not 0 or less than 0
        area[area <= 0] = 3.14 * 0.001**2
        return pd.DataFrame(
            {
                "kappa_area": params.shot_props.interpolation_method(
                    times, area / (np.pi * aminor**2), params.shot_props.times
                )
            }
        )


# --8<-- [end:kappa_area_request_example]
