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
        import os
        settings_file_path = os.path.expanduser("~/proj/ONW/ONE_SCOPE/disruption-py/disruption_py/config.toml")
        configs[tokamak] = Dynaconf(
            envvar_prefix="DISPY",
            # settings_file="disruption_py/config.toml",
            settings_file=settings_file_path,
            environments=True,
            default_env="default",
            env=tokamak,
            merge_enabled=True,
        )
    return configs[tokamak]
