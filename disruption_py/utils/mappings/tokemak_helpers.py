from disruption_py.databases import D3DDatabase, CModDatabase
from disruption_py.shots import D3DShot, CModShot
from disruption_py.utils.mappings.tokemak import Tokemak

DATABASE_HANDLERS = {Tokemak.D3D: D3DDatabase, Tokemak.CMOD: CModDatabase, Tokemak.EAST: None}
SHOT_CLASSES = {Tokemak.D3D: D3DShot, Tokemak.CMOD: CModShot, Tokemak.EAST: None}

def get_tokemak_from_shot_id(shot_id):
    if isinstance(shot_id, str):
        shot_len = len(shot_id)
    elif isinstance(shot_id, int):
        # math.log10 is faster and safer for large numbers but we assume shot_id is relatively small
        shot_len = len(str(shot_id))
    else:
        raise ValueError(
            f"shot_id must be a string or integer, not {type(shot_id)}")

    if shot_len == 6:
        return Tokemak.D3D
    elif shot_len == 10:
        return Tokemak.CMOD
    else:
        raise NotImplementedError(
            f"Unable to handle shot_id of length {shot_len}")

def get_database_for_shot_id(shot_id : int):
    tokemak = get_tokemak_from_shot_id(shot_id)
    return DATABASE_HANDLERS.get(tokemak, None)

def get_shot_class_for_shot_id(shot_id : int):
    tokemak = get_tokemak_from_shot_id(shot_id)
    return SHOT_CLASSES.get(tokemak, None)