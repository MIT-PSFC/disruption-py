#!/usr/bin/env python3

from typing import List

import pandas as pd

from disruption_py.shots.helpers.parameter_method_params import (
    ParameterMethodParams,
)
from disruption_py.shots.helpers.method_caching import parameter_method


def method_metadata_function(
    parent_object, parameter_method_params: ParameterMethodParams, **kwargs
) -> List[str]:
    """
    Parameters
    ----------
    parent_object : Any
        The object that contains the decorated method.
    parameter_method_params : ParameterMethodParams
        The ParameterMethodParams used when calling the decorated method.
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
    parent_object, parameter_method_params: ParameterMethodParams, **kwargs
) -> List[str]:
    # any properties of the `ShotProps` can be used to compute returned values
    if parameter_method_params.shot_props.shot_id > 10000000:
        return ["kappa_area", "q0"]
    else:
        return ["kappa_area"]


@parameter_method(columns=["kappa_area"])
def decorated_parameter_method(self, params: ParameterMethodParams) -> pd.DataFrame:
    pass


# --8<-- [end:decorator_functions_example]
