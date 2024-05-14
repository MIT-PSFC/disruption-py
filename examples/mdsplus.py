#!/usr/bin/env python3

"""
execute a simple fetch to test MDSplus connection.
"""

import os
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment, get_tokamak_handler

tokamak = get_tokamak_from_environment()
handler = get_tokamak_handler(tokamak)
if os.getenv("DIIID_TEST", False) or os.path.exists("/fusion/projects/disruption_warning"):
    shot = 161228
    shape = (196,)
else:
    shot = 1150805012
    shape = (2400,)
else:
    raise ValueError("Unspecified or unsupported tokamak.")

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
