#!/usr/bin/env python3

"""
This module provides functions to convert string attributes of an object
to corresponding enum values and to convert string values to enum values.
"""


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
        return value
    if value in enum_class:
        return enum_class(value)
    if value in enum_class.__members__:
        return enum_class[value]
    if should_raise:
        raise ValueError(f"Value '{value}' is not a valid member of {enum_class}.")
    return None
