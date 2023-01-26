"""
Check which shots have efit timebases that differ from database timebase and table by whether the shots are disrupted or not.
"""

from disruption_py.database import *

d3d = create_d3d_handler()
database_shots = d3d.get_disruption_table_shotlist()
for shot_id in database_shots:
    local_shot = d3d.get_shot(shot_id)
