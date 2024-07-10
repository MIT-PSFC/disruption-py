#!/usr/bin/env python3

"""
execute a simple fetch to test MDSplus connection.
"""

from disruption_py.machine.tokamak import Tokamak, get_tokamak_from_environment
from disruption_py.workflow import get_mdsplus_class

tokamak = get_tokamak_from_environment()

if tokamak is Tokamak.D3D:
    shot = 161228
    shape = (196,)
elif tokamak is Tokamak.CMOD:
    shot = 1150805012
    shape = (2400,)
else:
    raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

mds = get_mdsplus_class(tokamak).conn
print(f"Initialized MDSplus: {mds.hostspec}")

mds.openTree("EFIT01", shot)
print("#", shot)

node = r"\efit_aeqdsk:time"
print(">", node)

out = mds.get(node).data()
print("=", out.shape)
print(out)

assert out.shape == shape
