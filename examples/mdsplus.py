#!/usr/bin/env python3

"""
example module for MDSplus.
"""

from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.workflow import get_mdsplus_class


def main():
    """
    execute a simple fetch to test MDSplus connection.
    """

    tokamak = resolve_tokamak_from_environment()

    if tokamak is Tokamak.D3D:
        shot = 161228
        shape = (196,)
        efit = "efit01"
    elif tokamak is Tokamak.CMOD:
        shot = 1150805012
        shape = (62,)
        efit = "analysis"
    elif tokamak is Tokamak.EAST:
        shot = 55555
        shape = (102,)
        efit = "efit_east"
    else:
        raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

    mds = get_mdsplus_class(tokamak).conn
    print(f"Initialized MDSplus: {mds.hostspec}")

    mds.openTree(efit, shot)
    print("#", shot)

    node = r"dim_of(\efit_aeqdsk:li)"
    print(">", node)

    out = mds.get(node).data()
    print("=", out.shape)
    print(out)

    assert out.shape == shape


if __name__ == "__main__":
    main()
