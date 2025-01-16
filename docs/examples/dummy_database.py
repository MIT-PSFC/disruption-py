#!/usr/bin/env python3

"""
Example usage of `get_shots_data` demonstrating the use of a dummy sql database.
"""

from disruption_py.inout.sql import DummyDatabase
from disruption_py.workflow import get_shots_data

shot_data = get_shots_data(
    tokamak="d3d",
    shotlist_setting=[161228],
    database_initializer=DummyDatabase.initializer,
)
