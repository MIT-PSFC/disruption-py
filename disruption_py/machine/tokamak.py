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
    HBTEP = "hbtep"
    MAST = "mast"


def resolve_tokamak_from_environment(tokamak: Union[Tokamak, str] = None):
    """
    Method to resolve the tokamak:
    1. return if it's already a tokamak;
    2. read the argument, and overwrite the config;
    3. read the config;
    4. look for specific folders that will indicate presence on a given machine,
       and thus infer the tokamak from the cluster;
    5. return the only open-access tokamak as a fallback.
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

    tokamak_name = Tokamak.MAST  # default

    if tokamak:
        tokamak_name = map_string_to_enum(tokamak, Tokamak)
    # case 4
    if os.path.exists("/usr/local/mfe/disruptions"):
        tokamak_name = Tokamak.CMOD
    if os.path.exists("/fusion/projects/disruption_warning"):
        tokamak_name = Tokamak.D3D
    if os.path.exists("/project/disruption"):
        tokamak_name = Tokamak.EAST
    if os.path.exists("/opt/hbt/disruptions"):
        tokamak_name = Tokamak.HBTEP

    return tokamak_name
