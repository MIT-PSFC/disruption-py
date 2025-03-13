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
cache_dir = os.path.expanduser("~/.cache/disruption-py/")
os.makedirs(cache_dir, exist_ok=True)
cache_h5 = cache_dir + "my_workflow.h5"

tokamak = resolve_tokamak_from_environment()
shots = get_tokamak_test_shotlist(tokamak)

# Set log level to INFO to show the time spent on retrieving data
LOG_SETTINGS = "INFO"

initial_shots = shots[: len(shots) // 2]

# First run retrieves data for the first half of the shots
cached_shot_data = get_shots_data(
    shotlist_setting=initial_shots,
    log_settings=LOG_SETTINGS,
    output_setting=cache_h5,
)

# Second run to show little time is spent on retrieving data for the first half of the shots
rs = RetrievalSettings(
    cache_setting=cached_shot_data,
)
cached_shot_data = get_shots_data(
    retrieval_settings=rs,
    shotlist_setting=initial_shots,
    log_settings=LOG_SETTINGS,
)

# Final run on all the data, using the long-term saved cache data
retrieval_settings = RetrievalSettings(
    cache_setting=cache_h5,
)
shot_data = get_shots_data(
    tokamak=tokamak,
    shotlist_setting=shots,
    retrieval_settings=retrieval_settings,
    log_settings=LOG_SETTINGS,
)
