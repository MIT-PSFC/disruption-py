#!/usr/bin/env python3

"""Example usage of `get_shots_data` testing the connection to the SQL database."""

from disruption_py.workflow import get_database

cmod_database = get_database(tokamak="cmod")
shotlist = cmod_database.get_disruption_warning_shotlist()["shot"][0:8].tolist()
result = cmod_database.get_shots_data(shotlist, sql_table="disruption_warning")
