#!/usr/bin/env python3

"""
Example usage of `get_shots_data` when providing cache data from a previous run.
"""

import os

from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.utils.factory import get_tokamak_test_shotlist

# Create a long-lived cache directory, if it does not exist
# to demonstrate how to provide cache data from a previous run
cache_dir = os.path.expanduser("~/.cache/disruption-py")
os.makedirs(cache_dir, exist_ok=True)
cache_nc = os.path.join(cache_dir, "example_cache.nc")

tokamak = resolve_tokamak_from_environment()
shots = get_tokamak_test_shotlist(tokamak)

# Set log level to INFO to show the time spent on retrieving data
LOG_SETTINGS = "INFO"

initial_shots = shots[: len(shots) // 2]

# 1: retrieve and store data for some shots, takes some time
cached_shot_data = get_shots_data(
    shotlist_setting=initial_shots,
    output_setting=cache_nc,
    log_settings=LOG_SETTINGS,
)

# 2: use in-memory cache for the same shots, takes no time at all
cached_shot_data = get_shots_data(
    shotlist_setting=initial_shots,
    retrieval_settings=RetrievalSettings(cache_setting=cached_shot_data),
    log_settings=LOG_SETTINGS,
)

# 3: all the shots, priming cache from the long-term storage netcdf file
shot_data = get_shots_data(
    shotlist_setting=shots,
    tokamak=tokamak,
    retrieval_settings=RetrievalSettings(
        cache_setting=cache_nc,
    ),
    log_settings=LOG_SETTINGS,
)
