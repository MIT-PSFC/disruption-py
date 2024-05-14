#!/usr/bin/env python3

"""
execute a simple workflow to fetch EFIT parameters.
"""

from disruption_py.settings import ShotSettings, LogSettings
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment, get_tokamak_handler

tokamak = get_tokamak_from_environment()
handler = get_tokamak_handler(tokamak)
set_times_request = "disruption_warning"

if tokamak is Tokamak.D3D:
    shot_ids_request = [161228]
    run_methods = ["_get_efit_parameters"]
    shape = (247, 16)
elif tokamak is Tokamak.CMOD:
    shot_ids_request = [1150805012]
    run_methods = ["_get_EFIT_parameters"]
    shape = (62, 25)
else:
    raise ValueError(f"Tokamak {tokamak} not supported for this example")

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
