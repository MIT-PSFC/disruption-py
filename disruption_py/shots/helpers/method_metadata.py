#!/usr/bin/env python3

from dataclasses import Field, dataclass, fields
from typing import Any, Callable, List, Union

from disruption_py.settings.shot_data_request import ShotDataRequestParams
from disruption_py.utils.mappings.tokamak import Tokamak


@dataclass(frozen=True)
class MethodMetadata:
    name: str
    populate: bool

    cache_between_threads: bool
    used_trees: Union[List[str], Callable]
    contained_registered_methods: Union[List[str], Callable]
    tokamaks: Union[Tokamak, List[Tokamak]]
    columns: Union[List[str], Callable]
    tags: List[str]

    ALLOWED_UNRESOLVED = [
        "used_trees",
        "contained_registered_methods",
        "tokamaks",
    ]

    def __post_init__(self):
        if self.populate:
            object.__setattr__(self, "tags", self.tags or ["all"])
            object.__setattr__(self, "columns", self.columns or [])
        object.__setattr__(self, "used_trees", self.used_trees or [])
        object.__setattr__(self, "contained_registered_methods", [])
        object.__setattr__(self, "tokamaks", self.tokamaks or [])


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

        Some parameters provided to the register_method decorators can take method that are evaluated
        at runtime. `resolve_for` evaluates all of these methods and returns a new instance of `MethodMetadata`
        without function parameters.
        """
        new_method_metadata_params = {}
        bind_to = (getattr(bound_method, "__self__", None),)
        for field in fields(method_metadata):
            field_value = getattr(method_metadata, field.name)
            if field.name in method_metadata.ALLOWED_UNRESOLVED and callable(
                field_value
            ):
                new_val = (
                    field_value(params)
                    if bind_to is None
                    else field_value(bind_to, params)
                )
                new_method_metadata_params[field.name] = new_val
            else:
                new_method_metadata_params[field.name] = field_value

        return cls(bound_method=bound_method, **new_method_metadata_params)


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
    method : Callable
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
        raise ValueError(f"The method {method} was not decorated with register_method")
    return method_metadata
