#!/usr/bin/env python3

from typing import List

import pandas as pd

from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.physics_method.decorator import physics_method


def method_metadata_function(
    parent_object, physics_method_params: PhysicsMethodParams, **kwargs
) -> List[str]:
    """
    Parameters
    ----------
    parent_object : Any
        The object that contains the decorated method.
    physics_method_params : PhysicsMethodParams
        The PhysicsMethodParams used when calling the decorated method.
    kwargs : dict
        For future compatability.

    Returns
    -------
    List[str]
        Value of the expected type for the parameter (cuurently all parameters that allow functions are lists of strings)
    """
    pass


# --8<-- [start:decorator_functions_example]
def used_columns_by_shot_id(
    parent_object, physics_method_params: PhysicsMethodParams, **kwargs
) -> List[str]:
    # any properties of the `ShotProps` can be used to compute returned values
    if physics_method_params.shot_id > 10000000:
        return ["kappa_area", "q0"]
    else:
        return ["kappa_area"]


@physics_method(columns=["kappa_area"])
def decorated_physics_method(self, params: PhysicsMethodParams) -> pd.DataFrame:
    pass


# --8<-- [end:decorator_functions_example]
