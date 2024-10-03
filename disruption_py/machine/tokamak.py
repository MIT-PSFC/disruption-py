#!/usr/bin/env python3

"""Module for handling tokamak types and resolving tokamak configurations."""

import os
from enum import Enum
from typing import Union

from disruption_py.config import config
from disruption_py.core.utils.enums import map_string_to_enum


class Tokamak(Enum):
    """
    For documentation of supported tokamaks:
    # --8<-- [start:allowed_tokamak_types_snippet]
    Currently supported tokamak type strings are: `"cmod", "d3d"`
    # --8<-- [end:allowed_tokamak_types_snippet]
    """

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


def resolve_tokamak_from_environment(tokamak: Union[Tokamak, str] = None):
    """
    Method to resolve the tokamak:
    1. return if it's already a tokamak;
    2. read the argument, and overwrite the config;
    3. read the config;
    4. look for specific folders that will indicate presence on a given machine,
       and thus infer the tokamak from the cluster;
    5. raise exception.
    """
    if isinstance(tokamak, Tokamak):
        # case 1
        return tokamak
    if tokamak:
        # case 2
        config().tokamak = tokamak
    else:
        # case 3
        tokamak = config().get("tokamak")
    if tokamak:
        return map_string_to_enum(tokamak, Tokamak)
    # case 4
    if os.path.exists("/usr/local/mfe/disruptions"):
        return Tokamak.CMOD
    if os.path.exists("/fusion/projects/disruption_warning"):
        return Tokamak.D3D
    # case 5
    raise ValueError(
        "Tokamak is unspecified and could not be determined from the environment."
    )
