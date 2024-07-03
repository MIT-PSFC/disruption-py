#!/usr/bin/env python3

import functools
import threading
from typing import Callable, List

import numpy as np
import pandas as pd

from disruption_py.shots.helpers.parameter_method_params import ParameterMethodParams
from disruption_py.core.physics_method.params import ShotProps


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
