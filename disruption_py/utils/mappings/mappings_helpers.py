#!/usr/bin/env python3

from typing import Dict


def map_string_attributes_to_enum(obj, enum_translations: Dict):
    for field_name, enum_class in enum_translations.items():
        if hasattr(obj, field_name) and not isinstance(
            getattr(obj, field_name), enum_class
        ):
            try:
                enum_value = enum_class(getattr(obj, field_name))
                setattr(obj, field_name, enum_value)
            except ValueError:
                raise ValueError(
                    f"Invalid enum value for field '{field_name}': {getattr(obj, field_name)}"
                )


def map_string_to_enum(value, enum_class, should_raise=True):
    if isinstance(value, enum_class):
        enum_value = value
    else:
        try:
            enum_value = enum_class(value)
        except ValueError:
            if should_raise:
                raise ValueError(
                    f"Cannot convert value '{value}' to enum for '{enum_class}'"
                )
            else:
                return None
    return enum_value
