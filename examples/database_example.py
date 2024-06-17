#!/usr/bin/env python3

from disruption_py.database import ShotDatabase
from disruption_py.utils.constants import MDSPLUS_CONNECTION_STRING_CONSTANTS
from disruption_py.utils.mappings.tokamak import Tokamak


cmod_database = ShotDatabase.from_config(
    MDSPLUS_CONNECTION_STRING_CONSTANTS, tokamak=Tokamak.CMOD
)
shod_ids = cmod_database.get_disruption_warning_shotlist()["shot"][0:8].tolist()
result = cmod_database.get_shots_data(shod_ids, sql_table="disruption_warning")
