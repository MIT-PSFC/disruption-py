from disruption_py.shots.shot import Shot
from disruption_py.shots.d3d_shot import D3DShot, D3D_DISRUPTED_SHOT
from disruption_py.shots.cmod_shot import CModShot


def get_shot_id_type(shot_id):
    if isinstance(shot_id, str):
        shot_len = len(shot_id)
    elif isinstance(shot_id, int):
        # math.log10 is faster and safer for large numbers but we assume shot_id is relatively small
        shot_len = len(str(shot_id))
    else:
        raise ValueError(
            f"shot_id must be a string or integer, not {type(shot_id)}")

    if shot_len == 6:
        return D3DShot
    elif shot_len == 10:
        return CModShot
    else:
        raise NotImplementedError(
            f"Unable to handle shot_id of length {shot_len}")
