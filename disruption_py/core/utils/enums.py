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
