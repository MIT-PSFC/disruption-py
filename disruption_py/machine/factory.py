#!/usr/bin/env python3


from disruption_py.shots.cmod_shot_manager import CModShotManager
from disruption_py.shots.d3d_shot_manager import D3DShotManager
from disruption_py.machine.tokamak import Tokamak


def get_tokamak_shot_manager(tokamak: Tokamak):
    if tokamak == Tokamak.CMOD:
        return CModShotManager
    elif tokamak == Tokamak.D3D:
        return D3DShotManager
    else:
        raise ValueError("No shot manager for tokamak {}".format(tokamak))
