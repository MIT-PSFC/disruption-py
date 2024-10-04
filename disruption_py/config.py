#!/usr/bin/env python3

"""
Loads configuration settings using Dynaconf for a given tokamak.
"""

import os
from enum import Enum
from typing import Union

from dynaconf import Dynaconf

configs = {}


def config(tokamak: Union[Enum, str] = None):
    """
    Load and cache the configuration.

    Parameters
    ----------
    tokamak : Union[Enum, str], optional
        Tokamak name or Enum. Defaults to "default".

    Returns
    -------
    Dynaconf
        Configuration settings.
    """
    if tokamak is None:
        tokamak = "default"
    elif isinstance(tokamak, Enum):
        tokamak = tokamak.value

    if tokamak not in configs:
        configs[tokamak] = Dynaconf(
            envvar_prefix="DISPY",
            root_path=os.path.dirname(__file__),
            settings_file="config.toml",
            environments=True,
            default_env="default",
            env=tokamak,
            merge_enabled=True,
        )
    return configs[tokamak]
