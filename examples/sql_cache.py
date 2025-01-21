#!/usr/bin/env python3

"""
Example usage of `get_shots_data` when providing cache data."""

from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data

retrieval_settings = RetrievalSettings(
    # retrieve existing data from the disruption_warnings table use it to prepopulate queries
    cache_setting="sql",
)
shot_data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=[1150805012],
    retrieval_settings=retrieval_settings,
)
