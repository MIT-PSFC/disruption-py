#!/usr/bin/env python3

"""
Example usage of `get_shots_data` when providing cache data from a previous run.
"""

from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data

shots = [1150805012]
cached_shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=shots,
)

shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=shots,
    retrieval_settings=RetrievalSettings(
        cache_setting=cached_shot_data,
    ),
)
