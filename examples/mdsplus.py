#!/usr/bin/env python3

"""
execute a simple fetch to test MDSplus connection.
"""

from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment
from disruption_py.utils.eval.environment_constants import get_test_handler

tokamak = get_tokamak_from_environment()
if tokamak == Tokamak.D3D:
    shot = 161228
    shape = (196,)
elif tokamak == Tokamak.CMOD:
    shot = 1150805012
    shape = (2400,)
else:
    raise RuntimeError("Unspecified or unsupported tokamak.")

handler = get_test_handler(tokamak)
mds = handler.mds_connection.conn
print(f"Initialized MDSplus: {mds.hostspec}")

mds.openTree("EFIT01", shot)
print("#", shot)

node = r"\efit_aeqdsk:time"
print(">", node)

out = mds.get(node).data()
print("=", out.shape)
print(out)

assert out.shape == shape
