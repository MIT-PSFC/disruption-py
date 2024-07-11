#!/usr/bin/env python3

"""
execute a simple workflow to fetch EFIT parameters.
"""

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

tokamak = resolve_tokamak_from_environment()

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

retrieval_settings = RetrievalSettings(
    run_tags=[],
    run_methods=run_methods,
    efit_nickname_setting="analysis",
)

result = get_shots_data(
    tokamak=tokamak,
    shotlist_setting=shotlist,
    retrieval_settings=retrieval_settings,
    output_setting="dataframe",
    log_settings=LogSettings(console_log_level=0),
)

print(result)
assert result.shape == shape
