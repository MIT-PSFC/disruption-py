#!/usr/bin/env python3

import os
from enum import Enum
from logging import Logger

from disruption_py.core.utils.enums import map_string_to_enum

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


def get_tokamak_from_environment():
    if "DISPY_TOKAMAK" in os.environ:
        return Tokamak[os.environ["DISPY_TOKAMAK"]]
    if os.path.exists("/usr/local/mfe/disruptions"):
        return Tokamak.CMOD
    if os.path.exists("/fusion/projects/disruption_warning"):
        return Tokamak.D3D
    return None


def resolve_tokamak(tokamak: Tokamak, logger: Logger = None):
    if tokamak is None:
        tokamak = get_tokamak_from_environment()
        if logger:
            logger.info(f"No tokamak argument given. Detected tokamak: {tokamak.value}")
        return tokamak
    else:
        return map_string_to_enum(tokamak, Tokamak)
