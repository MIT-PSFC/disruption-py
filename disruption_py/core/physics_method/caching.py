#!/usr/bin/env python3

import functools
import threading
from typing import Callable, List

import numpy as np
import pandas as pd

from disruption_py.core.physics_method.params import PhysicsMethodParams


def cache_method(method: Callable) -> Callable:
    """Decorates a function as a cached method and instantiates its cache.

    Cached methods are functions that run expensive operations on data in the shot
    and may be reused. The cache is used to store the results of the parameter
    method so that it is only calculated once per shot for a given timebase.

    Parameters
    ----------
    method: Callable
        The method to be cached.
    """

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        physics_method_params: PhysicsMethodParams = (
            kwargs["params"] if "params" in kwargs else args[-1]
        )

        other_params = {k: v for k, v in kwargs.items() if k != "params"}

        cache_key = get_method_cache_key(
            method, physics_method_params.times, other_params
        )

        if cache_key in physics_method_params._cached_results:
            return physics_method_params._cached_results[cache_key]
        result = method(*args, **kwargs)
        physics_method_params._cached_results[cache_key] = result
        return result

    if isinstance(method, staticmethod):
        return staticmethod(wrapper)
    if isinstance(method, classmethod):
        return classmethod(wrapper)
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


def manually_cache(
    physics_method_params: PhysicsMethodParams,
    data: pd.DataFrame,
    method: Callable,
    method_name: str,
    method_columns: List[str],
) -> bool:
    if method_columns is None:
        return False
    if not hasattr(physics_method_params, "_cached_results"):
        physics_method_params._cached_results = {}
    missing_columns = set(col for col in method_columns if col not in data.columns)
    if len(missing_columns) == 0:
        cache_key = get_method_cache_key(method, data["time"].values)
        physics_method_params._cached_results[cache_key] = data[method_columns]
        physics_method_params.logger.debug(
            "[Shot %s]: Manually caching %s", physics_method_params.shot_id, method_name
        )
        return True
    physics_method_params.logger.debug(
        "[Shot %s]: Can not cache %s missing columns %s",
        physics_method_params.shot_id,
        method_name,
        missing_columns,
    )
    return False
