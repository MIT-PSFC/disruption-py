#!/usr/bin/env python3

"""
This module provides functions to convert string attributes of an object
to corresponding enum values and to convert string values to enum values.
"""

from typing import Dict


def map_string_attributes_to_enum(obj, enum_translations: Dict):
    """
    Map string attributes of an object to corresponding enum values.

    Parameters
    ----------
    obj : object
        The object whose attributes will be mapped to enum values.
    enum_translations : Dict[str, type]
        A dictionary mapping attribute names to their corresponding enum classes.
    """
    for field_name, enum_class in enum_translations.items():
        if hasattr(obj, field_name) and not isinstance(
            getattr(obj, field_name), enum_class
        ):
            try:
                enum_value = enum_class(getattr(obj, field_name))
                setattr(obj, field_name, enum_value)
            except ValueError as e:
                raise ValueError(
                    f"Invalid enum value for field '{field_name}': {getattr(obj, field_name)}"
                ) from e


def map_string_to_enum(value, enum_class, should_raise=True):
    """
    Convert a string value to the corresponding enum value.

    Parameters
    ----------
    value : str
        The string value to convert to an enum.
    enum_class : type
        The enum class to which the value should be converted.
    should_raise : bool, optional
        Whether to raise an exception if the conversion fails (default is True).

    Returns
    -------
    enum_class
        The corresponding enum value if conversion is successful, otherwise None
        if should_raise is False.
    """
    if isinstance(value, enum_class):
        enum_value = value
    else:
        try:
            enum_value = enum_class(value)
        except ValueError as e:
            if should_raise:
                raise ValueError(
                    f"Cannot convert value '{value}' to enum for '{enum_class}'"
                ) from e
            return None
    return enum_value
