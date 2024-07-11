from disruption_py.machine.d3d.basic import BasicD3DRequests
from disruption_py.machine.d3d.efit import D3DEfitRequests
from disruption_py.machine.cmod.basic import BasicCmodRequests, CModEfitRequests
from disruption_py.machine.tokamak import Tokamak


def get_physics_method_holders(tokamak: Tokamak):
    """Return a list of classes containing the built-in physics methods."""
    if tokamak is Tokamak.D3D:
        return [
            D3DEfitRequests(),
            BasicD3DRequests(),
        ]
    elif tokamak is Tokamak.CMOD:

        return [
            CModEfitRequests(),
            BasicCmodRequests(),
        ]
    else:
        raise ValueError(f"Invalid tokamak for physics methods {tokamak}")
