#!/usr/bin/env python3

import pandas as pd

from disruption_py.databases.database import ShotDatabase


class DummyDatabase(ShotDatabase):
    """
    A database class that does not require connecting to an SQL server but returns no data.

    Note: On CMod, disruption time data and any derrivative values will not be correct

    Examples
    --------
    >>> cmod_handler = CModHandler(database_initializer=DummyDatabase.default)
    >>> shot_data = cmod_handler.get_shots_data(shot_ids_request=[1150805012])
    <pd.DataFrame>
    """

    def __init__(self, **kwargs):
        pass

    @classmethod
    def default(cls, **kwargs):
        return cls()

    @property
    def conn(self, **kwargs):
        return DummyObject()

    def query(self, **kwargs):
        return pd.DataFrame()

    def get_shots_data(sefl, **kwargs):
        return pd.DataFrame()

    def get_disruption_time(self, **kwargs):
        return None

    def get_disruption_shotlist(self, **kwargs):
        return []

    def get_disruption_warning_shotlist(self, **kwargs):
        return []


class DummyObject:
    def __getattr__(self, name):
        # Return self for any attribute or method call
        return self

    def __call__(self, *args, **kwargs):
        # Return self for any method call
        return self
