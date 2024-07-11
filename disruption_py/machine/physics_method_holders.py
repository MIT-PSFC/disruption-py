#!/usr/bin/env python3

from disruption_py.machine.tokamak import Tokamak


def get_physics_method_holders(tokamak: Tokamak):
    """Return a list of classes containing the built-in physics methods."""
    if tokamak is Tokamak.D3D:
        from disruption_py.machine.d3d import METHOD_HOLDERS
        return METHOD_HOLDERS
    elif tokamak is Tokamak.CMOD:
        from disruption_py.machine.cmod import METHOD_HOLDERS
        return METHOD_HOLDERS
    else:
        raise ValueError(f"Invalid tokamak for physics methods {tokamak}")
