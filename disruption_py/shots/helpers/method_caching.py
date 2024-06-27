#!/usr/bin/env python3

from collections import defaultdict
import functools
import threading
from typing import Callable, List

import pandas as pd

from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.shots.helpers.method_metadata import (
    MethodMetadata,
)
from disruption_py.shots.shot_props import ShotProps

global_methods_registry: dict[str, list[MethodMetadata]] = defaultdict(list)
global_instances_registry = defaultdict(list)


def register_method(
    populate=True,
    cache=True,
    tokamak=None,
    tags=None,
    columns=None,
):
    """Registers a method to be run by DisruptionPy.

    All registered methods have their results cached, to avoid excess computation.
    If populate is True, the function should calculate disruption parameters and return a pandas DataFrame.
    A registered method with populate set to true is run if  designated by the run_methods, run_tags, or
    run_columns attributes of the `ShotSettings` class. A number of built-in methods are included
    through the shots.parameter_methods.built_in file. Users can also create there own methods
    using the register_method , and passing either the method itself, or an object containing the method
    as and attribute in `ShotSettings`.

    A common pattern for parameter methods is first retrieving data from MDSplus using the `TreeManager` and
    then using that retrieved data to compute data to return.

    Parameters
    ----------
    populate : bool
        Whether this method should be run when calling `get_shots_data`. If True the method must return a dictionary of
        column name to data. Default is True.
    tags : list[str]
        The list of tags to help identify this method. Tags can be used when calling `get_shots_data` to
        have disruption_py use multiple functions at the same time. Default is ["all"].
    columns : Union[list[str], Callable]
        The columns that are in the dataframe returned by the method. Alternately, can pass a method that
        returns the names of used trees at runtime. See `method_metadata_function` for more details about using functions.
        These columns are also used to determine which methods to run when calling `get_shots_data` with `run_columns`.
        Default value is an empty list implying that no columns are returned by the function.
    cache: bool
        Whether to cache the result of the method. Default is True.
    tokamak : Union[Tokamak, List[Tokamak]]
        A list of Tokamak objects that represent which tokamaks this method may be used for. Default value of None allows the method
        to be run for any tokamak.
    """

    def outer_wrapper(method):
        @functools.wraps(method)
        def wrapper(*args, **kwargs):
            if cache:
                return cache_or_compute(method, args, kwargs)
            else:
                return method(*args, **kwargs)

        method_metadata = MethodMetadata(
            name=method.__name__,
            cache=cache,
            populate=populate,
            tokamaks=tokamak,
            columns=columns,
            tags=tags,
        )

        wrapper.method_metadata = method_metadata

        # Register the method in the global registry
        if isinstance(method, staticmethod):
            return staticmethod(wrapper)
        elif isinstance(method, classmethod):
            return classmethod(wrapper)
        else:
            return wrapper

    return outer_wrapper


def get_method_cache_key(method_name, times):
    current_thread_id = threading.get_ident()
    return method_name + str(len(times)) + str(current_thread_id)


def cache_or_compute(
    func: Callable,
    args: tuple,
    kwargs: dict,
):
    shot_data_request_params: ShotDataRequestParams = (
        kwargs["params"] if "params" in kwargs else args[-1]
    )
    shot_props = shot_data_request_params.shot_props

    other_params = {k: v for k, v in kwargs.items() if k != "params"}

    cache_key = get_method_cache_key(func.__name__, shot_props.times) + str(
        other_params
    )
    if cache_key in shot_props._cached_results:
        return shot_props._cached_results[cache_key]
    else:
        result = func(*args, **kwargs)
        shot_props._cached_results[cache_key] = result
        return result


def manually_cache(
    shot_props: ShotProps, data: pd.DataFrame, method_name, method_columns: List[str]
) -> bool:
    if method_columns is None:
        return False
    if not hasattr(shot_props, "_cached_results"):
        shot_props._cached_results = {}
    missing_columns = set(col for col in method_columns if col not in data.columns)
    if len(missing_columns) == 0:
        cache_key = get_method_cache_key(method_name, data["time"])
        shot_props._cached_results[cache_key] = data[method_columns]
        shot_props.logger.debug(
            f"[Shot {shot_props.shot_id}]:Manually caching {method_name}"
        )
        return True
    else:
        shot_props.logger.debug(
            f"[Shot {shot_props.shot_id}]:Can not cache {method_name} missing columns {missing_columns}"
        )
        return False
