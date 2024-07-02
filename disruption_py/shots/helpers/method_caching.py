#!/usr/bin/env python3

import functools
import threading
from typing import Callable, List, Union

import numpy as np
import pandas as pd

from disruption_py.shots.helpers.parameter_method_params import ParameterMethodParams
from disruption_py.shots.helpers.method_metadata import (
    MethodMetadata,
)
from disruption_py.shots.shot_props import ShotProps
from disruption_py.machine.tokamak import Tokamak


def parameter_method(
    cache: bool = True,
    tokamak: Union[Tokamak, List[Tokamak]] = None,
    tags: List[str] = None,
    columns: Union[List[str], Callable] = None,
) -> Callable:
    """
    Decorator to signify a method to be run by DisruptionPy.

    The decorated method calculates disruption paramaters and returns a pandas DataFrame. All decorated methods must take
    the single argument params of type `ParameterMethodParams`. The decorated method will be run if designated by the
    `run_methods`, `run_tags`, or `run_columns` attributes of the `ShotSettings` class, and if included inside of the
    `parameter_methods` argument of the `shot_settings` or in the built-in method list. If run the result of the decorated
    method will be output to the `output_type_request`.

    A common pattern for parameterized methods is first retrieving data from MDSplus using the `TreeManager` and
    then using that retrieved data to compute data to return.

    Parameters
    ----------
    cache : bool, optional
        Whether to cache the result of the method, by default True.
    tokamak : Union['Tokamak', List['Tokamak']], optional
        A list of Tokamak objects that represent which tokamaks this method may be used for, by default None
        allows the method to be run for any tokamak.
    tags : List[str], optional
        The list of tags to help identify this method. Tags can be used when calling `get_shots_data` to
        have disruption_py use multiple functions at the same time. Default is ["all"].
    columns : Union[List[str], Callable], optional
        The columns that are in the DataFrame returned by the method. Alternately, can pass a method that
        returns the names of used trees at runtime. Default value is an empty list implying that no columns
        are returned by the function.

    Returns
    -------
    Callable
        A decorated method with caching and metadata attributes.
    """

    def outer_wrapper(method: Callable) -> Callable:
        if cache:
            wrapper = cache_method(method)
        else:
            wrapper = method

        method_metadata = MethodMetadata(
            name=method.__name__,
            cache=cache,
            populate=True,
            tokamaks=tokamak,
            columns=columns,
            tags=tags,
        )

        wrapper.method_metadata = method_metadata
        return wrapper

    return outer_wrapper


def cache_method(method: Callable) -> Callable:
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        return cache_or_compute(method, args, kwargs)

    # parameter the method in the global registry
    if isinstance(method, staticmethod):
        return staticmethod(wrapper)
    elif isinstance(method, classmethod):
        return classmethod(wrapper)
    else:
        return wrapper


def get_method_cache_key(
    method: Callable, times: np.ndarray, other_params: dict = None
):
    current_thread_id = threading.get_ident()
    hashable_other_params = frozenset((other_params or {}).items())
    return (
        current_thread_id,
        method,
        times[0],
        times[-1],
        len(times),
        hashable_other_params,
    )


def cache_or_compute(
    method: Callable,
    args: tuple,
    kwargs: dict,
):
    parameter_method_params: ParameterMethodParams = (
        kwargs["params"] if "params" in kwargs else args[-1]
    )
    shot_props = parameter_method_params.shot_props

    other_params = {k: v for k, v in kwargs.items() if k != "params"}

    cache_key = get_method_cache_key(method, shot_props.times, other_params)

    if cache_key in shot_props._cached_results:
        return shot_props._cached_results[cache_key]
    else:
        result = method(*args, **kwargs)
        shot_props._cached_results[cache_key] = result
        return result


def manually_cache(
    shot_props: ShotProps,
    data: pd.DataFrame,
    method: Callable,
    method_name: str,
    method_columns: List[str],
) -> bool:
    if method_columns is None:
        return False
    if not hasattr(shot_props, "_cached_results"):
        shot_props._cached_results = {}
    missing_columns = set(col for col in method_columns if col not in data.columns)
    if len(missing_columns) == 0:
        cache_key = get_method_cache_key(method, data["time"].values)
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
