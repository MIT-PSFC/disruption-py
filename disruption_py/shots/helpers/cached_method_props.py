#!/usr/bin/env python3

from dataclasses import Field, dataclass
from typing import Any, Callable, List, Union

from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass(frozen=True)
class MethodMetadata:
    name: str
    populate: bool

    cache_between_threads: bool
    used_trees: Union[List[str], Callable]
    contained_cached_methods: Union[List[str], Callable]
    tokamaks: Union[Tokamak, List[Tokamak]]
    columns: Union[List[str], Callable] = None
    tags: List[str] = None

    ALLOWED_UNRESOLVED = [
        "used_trees",
        "contained_cached_methods",
        "tokamaks",
    ]

    def __post_init__(self):
        if self.populate:
            self.tags = self.tags or ["all"]
            self.columns = self.columns or []


@dataclass(frozen=True)
class BoundMethodMetadata(MethodMetadata):
    bound_method: Callable

    @classmethod
    def bind(
        cls,
        method_metadata: MethodMetadata,
        bound_method: Callable,
        params: ShotDataRequestParams,
    ):
        """
        Evaluate arguments to decorators to usable values.

        Some parameters provided to the cached_method and parameter_cached_method decorators can take method that are evaluated
        at runtime. `resolve_for` evaluates all of these methods and returns a new instance of `MethodMetadata`
        without function parameters.
        """
        new_method_metadata_params = {}
        bind_to = (getattr(bound_method, "__self__", None),)
        for field_name in method_metadata.ALLOWED_UNRESOLVED:
            field_value = getattr(method_metadata, field_name)
            if callable(field_value):
                new_val = (
                    field_value(params)
                    if bind_to is None
                    else field_value(bind_to, params)
                )
                new_method_metadata_params[field_name] = new_val
            else:
                new_method_metadata_params[field_name] = field_value

        return cls(bound_method=bound_method, **new_method_metadata_params)


@dataclass
class CachedMethodProps:
    name: str
    method: Callable

    # All functions have been evaluated
    computed_method_metadata: MethodMetadata

    def get_param_value(self, field_name: str, default_value: Any = None) -> Any:
        return (
            getattr(self.computed_method_metadata, field_name, default_value)
            or default_value
        )


# Utility methods for decorated methods


def is_registered_method(method: Callable) -> bool:
    """Returns whether the method is decorated with `register_method` decorator

    Parameters
    ----------
    method : Callable
        The method to check if decorated.

    Returns
    -------
    bool
        Whether the passed method is decorated.
    """
    return hasattr(method, "method_metadata")


def get_method_metadata(method: Callable, should_throw: bool = False) -> MethodMetadata:
    """Get method metadata for method

    Parameters
    ----------
    cached_method : Callable
        The method decorated with the `register_method` decorator
    should_throw : bool
        Throw an error if the method was not decorated with the `register_method` decorator

    Returns
    -------
    MethodMetadata
        The `MethodMetadata` object holding parameters for the cached method
    """
    method_metadata = getattr(method, "method_metadata", None)
    if should_throw and method_metadata is None:
        raise ValueError(
            f"The method {method} was not decorated with cached_method or parameter_cached_method"
        )
    return method_metadata
