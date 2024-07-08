#!/usr/bin/env python3

from enum import Enum

"""
For documentation of supported tokamaks:
# --8<-- [start:allowed_tokamak_types_snippet]
Currently supported tokamak type strings are: `"cmod"`
# --8<-- [end:allowed_tokamak_types_snippet]
"""


class Tokamak(Enum):
    D3D = "d3d"
    CMOD = "cmod"
    EAST = "east"


def is_tokamak_indexed(check_dict: dict):
    """
    Check if a dictionary is indexed by tokamak.
    """
    for option in Tokamak:
        if option.value in check_dict:
            return True

    return False
