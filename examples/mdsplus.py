#!/usr/bin/env python3

"""
execute a simple fetch to test MDSplus connection.
"""

from disruption_py.mdsplus_integration.mds_connection import ProcessMDSConnection
from disruption_py.utils.constants import MDSPLUS_CONNECTION_STRING_CONSTANTS
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak import (
    get_tokamak_from_environment,
)

tokamak = get_tokamak_from_environment()

if tokamak is Tokamak.D3D:
    shot = 161228
    shape = (196,)
elif tokamak is Tokamak.CMOD:
    shot = 1150805012
    shape = (2400,)
else:
    raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

mds = ProcessMDSConnection.from_config(tokamak=tokamak).conn
print(f"Initialized MDSplus: {mds.hostspec}")

mds.openTree("EFIT01", shot)
print("#", shot)

node = r"\efit_aeqdsk:time"
print(">", node)

out = mds.get(node).data()
print("=", out.shape)
print(out)

assert out.shape == shape
