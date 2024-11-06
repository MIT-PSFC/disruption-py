#!/usr/bin/env python3

"""Module for retrieving physics method classes based on the specified tokamak."""

from disruption_py.machine.tokamak import Tokamak


def get_method_holders(tokamak: Tokamak):
    """
    Return a list of classes containing the built-in physics methods.
    """
    # such an import pattern lets us avoid dynamic imports or code introspection
    # pylint: disable=import-outside-toplevel
    from disruption_py.machine.generic import METHOD_HOLDERS as GENERIC

    if tokamak is Tokamak.D3D:
        from disruption_py.machine.d3d import METHOD_HOLDERS

        return GENERIC + METHOD_HOLDERS
    if tokamak is Tokamak.CMOD:
        from disruption_py.machine.cmod import METHOD_HOLDERS

        return GENERIC + METHOD_HOLDERS
    raise ValueError(f"Invalid tokamak for physics methods {tokamak}")
