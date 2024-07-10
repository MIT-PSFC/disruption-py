#!/usr/bin/env python3

import os
from enum import Enum

from disruption_py.config import config
from disruption_py.core.utils.enums import map_string_to_enum

"""
For documentation of supported tokamaks:
# --8<-- [start:allowed_tokamak_types_snippet]
Currently supported tokamak type strings are: `"cmod", "d3d"`
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


def resolve_tokamak_from_environment():
    """
    Method to resolve the tokamak.
    First, try reading from the configuration, then try looking for specific folders that will
    indicate presence on a given machine, and thus infer the tokamak from the cluster.
    """
    tokamak = config().get("tokamak")
    if tokamak:
        return map_string_to_enum(tokamak, Tokamak)
    if os.path.exists("/usr/local/mfe/disruptions"):
        return Tokamak.CMOD
    if os.path.exists("/fusion/projects/disruption_warning"):
        return Tokamak.D3D
    raise ValueError(
        "Tokamak is unspecified and could not be determined from the environment."
    )
