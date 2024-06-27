#!/usr/bin/env python3

from typing import List

import pandas as pd

from disruption_py.settings.shot_data_request import (
    ShotDataRequestParams,
)
from disruption_py.shots.helpers.method_caching import register_method


def method_metadata_function(
    parent_object, shot_data_request_params: ShotDataRequestParams, **kwargs
) -> List[str]:
    """
    Parameters
    ----------
    parent_object : Any
        The object that contains the decorated method.
    shot_data_request_params : ShotDataRequestParams
        The ShotDataRequestParams used when calling the decorated method.
    kwargs : dict
        For future compatability.

    Returns
    -------
    List[str]
        Value of the expected type for the parameter (cuurently all parameters that allow functions are lists of strings)
    """
    pass


# --8<-- [start:decorator_functions_example]
def used_trees_by_shot_id(
    parent_object, shot_data_request_params: ShotDataRequestParams, **kwargs
) -> List[str]:
    # any properties of the `ShotProps` can be used to compute returned values
    if shot_data_request_params.shot_props.shot_id > 10000000:
        return ["tree_1", "tree_2"]
    else:
        return ["tree_1", "tree_3"]


@register_method(used_trees=["tree_1", "tree_2"])
def decorated_shot_data_method(self, params: ShotDataRequestParams) -> pd.DataFrame:
    pass


# --8<-- [end:decorator_functions_example]
