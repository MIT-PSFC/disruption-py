#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Any, Callable, List, Union

from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass(frozen=True)
class CachedMethodParams:
    cache_between_threads: bool
    used_trees: Union[List[str], Callable]
    contained_cached_methods: Union[List[str], Callable]
    tokamaks: Union[Tokamak, List[Tokamak]]


@dataclass(frozen=True)
class ParameterCachedMethodParams(CachedMethodParams):
    columns: Union[List[str], Callable]
    tags: List[str]

    def from_cached_method_params(
        cached_method_params: CachedMethodParams, columns, tags
    ) -> "ParameterCachedMethodParams":
        return ParameterCachedMethodParams(
            cache_between_threads=cached_method_params.cache_between_threads,
            used_trees=cached_method_params.used_trees,
            contained_cached_methods=cached_method_params.contained_cached_methods,
            tokamaks=cached_method_params.tokamaks,
            columns=columns,
            tags=tags,
        )


@dataclass
class CachedMethodProps:
    name: str
    method: Callable

    # All functions have been evaluated
    computed_cached_method_params: CachedMethodParams

    def get_param_value(self, field_name: str, default_value: Any = None) -> Any:
        return (
            getattr(self.computed_cached_method_params, field_name, default_value)
            or default_value
        )


# Utility methods for decorated methods


def is_cached_method(cached_method: Callable) -> bool:
    """Returns whether the method is decorated with `cached_method` or `parameter_cached_method` decorators

    Parameters
    ----------
    cached_method : Callable
        The method to check if decorated.

    Returns
    -------
    bool
        Whether the passed method is decorated.
    """
    return hasattr(cached_method, "cached_method_params")


def get_cached_method_params(
    cached_method: Callable, should_throw: bool = False
) -> CachedMethodParams:
    """Get cached method params for cached method

    Parameters
    ----------
    cached_method : Callable
        The method decorated with `cached_method` or `parameter_cached_method` decorators
    should_throw : bool
        Throw an error if the method was not decorated with `cached_method` or `parameter_cached_method` decorators

    Returns
    -------
    CachedMethodParams
        The `CachedMethodParams` object holding parameters for the cached method
    """
    cached_method_params = getattr(cached_method, "cached_method_params", None)
    if should_throw and cached_method_params is None:
        raise ValueError(
            f"The method {cached_method} was not decorated with cached_method or parameter_cached_method"
        )
    return cached_method_params
