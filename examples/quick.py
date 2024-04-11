#!/usr/bin/env python3

"""
quick test: fetch efit parameters from established shot, and check array shape.
"""

from disruption_py.handlers import CModHandler
from disruption_py.settings import ShotSettings, LogSettings

handler = CModHandler()

shot_settings = ShotSettings(
    log_settings=LogSettings(console_log_level=0),
    run_tags=[],
    run_methods=["_get_EFIT_parameters"],
    use_hsds=True,
    use_mongo=False,
    cache_miss_enable=False, 
)

result = handler.get_shots_data(
    shot_ids_request=[1150805012],
    shot_settings=shot_settings,
    output_type_request="dataframe",
)

print(result)

assert result.shape == (62, 25)
