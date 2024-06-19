#!/usr/bin/env python3

from disruption_py.main import get_database
from disruption_py.utils.mappings.tokamak import Tokamak


cmod_database = get_database(tokamak="cmod")
shod_ids = cmod_database.get_disruption_warning_shotlist()["shot"][0:8].tolist()
result = cmod_database.get_shots_data(shod_ids, sql_table="disruption_warning")
