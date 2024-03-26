#!/usr/bin/env python3

"""
execute a simple fetch to test MDSplus connection.
"""

import os
from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.handlers.d3d_handler import D3DHandler

if os.path.exists("/fusion/projects/disruption_warning"):
    handler = D3DHandler()
    shot = 161228
else:
    handler = CModHandler()
    shot = 1150805012
mds = handler.mds_connection.conn
print(f"Initialized MDSplus: {mds.hostspec}")

mds.openTree("EFIT01", shot)
print("#", shot)

node = r"\efit_aeqdsk:time"
print(">", node)

out = mds.get(node).data()
print("=", out.shape)
print(out)
