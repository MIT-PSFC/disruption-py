#!/usr/bin/env python3

from disruption_py.machine.tokamak import Tokamak


def built_in_method_factory(tokamak: Tokamak):
    if tokamak is Tokamak.D3D:
        from disruption_py.machine.d3d.basic import BasicD3DRequests
        from disruption_py.machine.d3d.efit import D3DEfitRequests

        return [
            D3DEfitRequests(),
            BasicD3DRequests(),
        ]
    elif tokamak is Tokamak.CMOD:
        from disruption_py.machine.cmod.basic import BasicCmodRequests, CModEfitRequests

        return [
            CModEfitRequests(),
            BasicCmodRequests(),
        ]
    else:
        raise ValueError(f"Invalid tokamak for built-ins {tokamak}")
