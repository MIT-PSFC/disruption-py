#!/usr/bin/env python3

"""
Module for defining metadata classes for physics methods.
"""

from dataclasses import dataclass, fields
from typing import Callable, List, Union

from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.tokamak import Tokamak


@dataclass(frozen=True)
class MethodMetadata:
    """
    Holder for the arguments to the decorator.
    """

    name: str

    cache: bool
    tokamaks: Union[Tokamak, List[Tokamak]]
    columns: Union[List[str], Callable]
    tags: List[str]

    ALLOWED_UNRESOLVED = [
        "columns",
        "tokamaks",
    ]

    def __post_init__(self):
        object.__setattr__(self, "tags", self.tags or ["all"])
        object.__setattr__(self, "columns", self.columns or [])


@dataclass(frozen=True)
class BoundMethodMetadata(MethodMetadata):
    """
    Metadata for a bound method, extending `MethodMetadata`.

    Attributes
    ----------
    bound_method : Callable
        The method that is bound to this metadata.
    """

    bound_method: Callable

    @classmethod
    def bind(
        cls,
        method_metadata: MethodMetadata,
        bound_method: Callable,
        physics_method_params: PhysicsMethodParams,
    ):
        """
        Bind a method to its metadata and resolve any callable parameters.

        Parameters
        ----------
        method_metadata : MethodMetadata
            Metadata instance containing the method's unresolved parameters.
        bound_method : Callable
            The callable method to be bound.
        physics_method_params : PhysicsMethodParams
            Parameters required for resolving the method.

        Returns
        -------
        BoundMethodMetadata
            A new instance of `BoundMethodMetadata` with resolved parameters.
        """
        new_method_metadata_params = {}
        bind_to = (getattr(bound_method, "__self__", None),)
        for field in fields(method_metadata):
            field_value = getattr(method_metadata, field.name)
            if field.name in method_metadata.ALLOWED_UNRESOLVED and callable(
                field_value
            ):
                new_val = (
                    field_value(physics_method_params)
                    if bind_to is None
                    else field_value(bind_to, physics_method_params)
                )
                new_method_metadata_params[field.name] = new_val
            else:
                new_method_metadata_params[field.name] = field_value

        return cls(bound_method=bound_method, **new_method_metadata_params)


# Utility methods for decorated methods


def is_parametered_method(method: Callable) -> bool:
    """
    Returns whether the method is decorated with `physics_method` decorator

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
    """
    Get method metadata for method

    Parameters
    ----------
    method : Callable
        The method decorated with the `physics_method` decorator
    should_throw : bool
        Throw an error if the method was not decorated with the `physics_method` decorator

    Returns
    -------
    MethodMetadata
        The `MethodMetadata` object holding parameters for the cached method
    """
    method_metadata = getattr(method, "method_metadata", None)
    if should_throw and method_metadata is None:
        raise ValueError(f"The method {method} was not decorated with physics_method")
    return method_metadata
