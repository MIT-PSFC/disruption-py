#!/usr/bin/env python3

from enum import Enum
from typing import Union

from dynaconf import Dynaconf

configs = {}


def config(tokamak: Union[Enum, str] = None):
    if tokamak is None:
        tokamak = "default"
    elif isinstance(tokamak, Enum):
        tokamak = tokamak.value

    if tokamak not in configs:
        configs[tokamak] = Dynaconf(
            envvar_prefix="DISPY",
            settings_file="disruption_py/config.toml",
            environments=True,
            default_env="default",
            env=tokamak,
            merge_enabled=True,
        )
    return configs[tokamak]
