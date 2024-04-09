import os
from disruption_py.databases import D3DDatabase, CModDatabase
from disruption_py.utils.mappings.tokamak import Tokamak

DATABASE_HANDLERS = {Tokamak.D3D: D3DDatabase, Tokamak.CMOD: CModDatabase, Tokamak.EAST: None}

def get_tokamak_from_shot_id(shot_id):
    if isinstance(shot_id, str):
        shot_len = len(shot_id)
    elif isinstance(shot_id, int):
        # math.log10 is faster and safer for large numbers but we assume shot_id is relatively small
        shot_len = len(str(shot_id))
    else:
        raise ValueError(
            f"shot_id must be a string or integer, not {type(shot_id)}")

    if shot_len == 6:
        return Tokamak.D3D
    elif shot_len == 10:
        return Tokamak.CMOD
    else:
        raise NotImplementedError(
            f"Unable to handle shot_id of length {shot_len}")

def get_tokamak_from_environment():
    if os.path.exists("/usr/local/mfe/disruptions"):
        return Tokamak.CMOD
    if os.path.exists("/fusion/projects/disruption_warning"):
        return Tokamak.D3D
    if "DISPY_TOKAMAK" in os.environ:
        return Tokamak[os.environ["DISPY_TOKAMAK"]]
    return None

def get_database_for_shot_id(shot_id : int):
    tokamak = get_tokamak_from_shot_id(shot_id)
    return DATABASE_HANDLERS.get(tokamak, None)
