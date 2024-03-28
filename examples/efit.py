#!/usr/bin/env python3

"""
execute a simple workflow to fetch EFIT parameters.
"""

import os
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.handlers.d3d_handler import D3DHandler
from disruption_py.settings import ShotSettings, LogSettings


if os.path.exists("/fusion/projects/disruption_warning"):
    handler = D3DHandler()
    shot_ids_request = [161228]
    set_times_request = "disruption"
    run_methods = ["_get_efit_parameters"]
    shape = (247, 16)
else:
    handler = CModHandler()
    shot_ids_request = [1150805012]
    set_times_request = "efit"
    run_methods = ["_get_EFIT_parameters"]
    shape = (62, 25)
print(f"Initialized handler: {handler.get_tokamak().value}")

shot_settings = ShotSettings(
    set_times_request=set_times_request,
    log_settings=LogSettings(console_log_level=0),
    run_tags=[],
    run_methods=run_methods,
)

result = handler.get_shots_data(
    shot_ids_request=shot_ids_request,
    shot_settings=shot_settings,
    output_type_request="dataframe",
)

print(result)

assert result.shape == shape
