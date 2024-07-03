#!/usr/bin/env python3

"""
execute a simple workflow to fetch EFIT parameters.
"""

from disruption_py.workflow import get_shots_data
from disruption_py.settings import LogSettings, Settings
from disruption_py.machine.tokamak import Tokamak
from disruption_py.machine.tokamak import (
    get_tokamak_from_environment,
)

tokamak = get_tokamak_from_environment()

if tokamak is Tokamak.D3D:
    shotlist = [161228]
    run_methods = ["_get_efit_parameters"]
    shape = (247, 16)
elif tokamak is Tokamak.CMOD:
    shotlist = [1150805012]
    run_methods = ["_get_EFIT_parameters"]
    shape = (62, 25)
else:
    raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

print(f"Initialized for tokamak: {tokamak.value}")

shot_settings = Settings(
    run_tags=[],
    run_methods=run_methods,
)

result = get_shots_data(
    tokamak=tokamak,
    shotlist_setting=shotlist,
    shot_settings=shot_settings,
    output_setting="dataframe",
    log_settings=LogSettings(console_log_level=0),
)

print(result)

assert result.shape == shape
