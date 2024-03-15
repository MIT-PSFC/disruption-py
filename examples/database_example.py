from disruption_py.handlers.cmod_handler import CModHandler

cmod_database = CModHandler().database
shot_ids = cmod_database.get_disruption_warning_shotlist()['shot'][0:8].tolist()
result = cmod_database.get_shots_data(shot_ids, sql_table="disruption_warning")