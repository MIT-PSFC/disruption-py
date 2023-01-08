"""
While trying to open efit trees for shots 191914 and 191786, I got the following error:
File "/fusion/projects/disruption_warning/disruption-warning-db-workflow/disruption_py/shots/d3d_shot.py", line 52, in get_efit_parameters
    self.conn.openTree(self.efit_tree_name, self._shot_id)
  File "/fusion/usc/src/mdsplus/mdsplus-stable_release-6-1-84/mdsobjects/python/build/lib/MDSplus/connection.py", line 211, in openTree
    status=self.get("TreeOpen($,$)",tree,shot)
  File "/fusion/usc/src/mdsplus/mdsplus-stable_release-6-1-84/mdsobjects/python/build/lib/MDSplus/connection.py", line 262, in get
    return self.__getAnswer__()
  File "/fusion/usc/src/mdsplus/mdsplus-stable_release-6-1-84/mdsobjects/python/build/lib/MDSplus/connection.py", line 132, in __getAnswer__
    raise MdsException(str(ans))
MDSplus._mdsshr.MdsException: %TDI-E-INVDTYDSC, Storage data type is not valid

This script generates a CSV containing all the shots that have this error for 'efit01'.
"""
from disruption_py.database import create_d3d_handler
import MDSplus


def find_shots_with_bad_trees(table='disruption_warning', tree_to_check=['efit01']):
    """
    This function connects to the MDSplus server, and loops over all shots in the
    specified DIII-D table, attempting to open the 'trees_to_check' trees for each shot.
    It then returns a list of shots for which at least one of the trees could not be opened, along with the
    number of shots that were attempted.

    Inputs:
        trees_to_check: the name of the trees to check. Default is just 'efit01'.

    Outputs:
        bad_shots: list of tuples containing the shot number and the error
            message from MDSplus for each shot that could not be opened.
        num_shots: the number of shots that were attempted
    """
    d3d = create_d3d_handler()
    bad_shots = []
    if table == 'disruption_warning':
        all_shots = d3d.get_disruption_table_shotlist()
    elif table == 'disruptions':
        all_shots = d3d.get_disruption_shotlist()
    else:
        raise ValueError(f"{table} table is not supported by the d3d handler")
    num_shots = len(all_shots)
    conn = MDSplus.Connection('atlas.gat.com')
    for shot in all_shots['shot']:
        for tree in tree_to_check:
            try:
                conn.openTree(tree, shot)
            except MDSplus._mdsshr.MdsException as e:
                bad_shots.append((shot, e))
                break
    return bad_shots, num_shots


if __name__ == '__main__':
    bad_shots, num_shots = find_shots_with_bad_trees('disruptions')
    with open('bad_shots_all.csv', 'w') as f:
        f.write('shot,error')
        for shot, error in bad_shots:
            f.write('{},{}\n'.format(shot, error))
    print(
        f"Found a total of {len(bad_shots)} bad shots(out of {num_shots}).")
