from MDSplus import *

DEFAULT_SHOT_COLUMNS = ['time','shot','time_until_disrupt','ip']

class Shot:
    """
    Base shot class to represent a single shot of a fusion experiment.

    Attributes
    ---
    metadata->dict: 
        labels: {} /* Dictionary where key is the label name(e.g. "disruption period")
                    and value is a list of starting and ending time tuples(e.g.[(10,100)]*/
        commit_hash: // Commit hash of code version that was used to generate shot data 
        timestep: // Time delta. Every array must conform to this time step
        duration: // Length in milliseconds of shot
        description: // Optional text summary to add to shot
    data->dict:
        array1:
        array2:
        ...
    """
    def __init__(self, mdsplus_name,shot_id,data_columns = DEFAULT_SHOT_COLUMNS):
        t = Tree(mdsplus_name,shot_id)


class D3DShot(Shot):
    pass

class EASTShot(Shot):
    pass

class CMODShot(Shot):
    pass 
    