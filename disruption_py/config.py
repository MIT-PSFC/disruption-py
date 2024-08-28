#!/usr/bin/env python3
import os
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
        #settings_file_path = os.path.expanduser("~/proj/ONW/ONE_SCOPE/disruption-py/disruption_py/config.toml")
        configs[tokamak] = Dynaconf(
            envvar_prefix="DISPY",
            #settings_file=settings_file_path, #TODO(ZanderKeith): Submit an issue about this relative path not working
            environments=True,
            default_env="default",
            env=tokamak,
            merge_enabled=True,
        )
    return configs[tokamak]
