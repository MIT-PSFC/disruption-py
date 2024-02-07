import sys
sys.path.append("/home/joshlor/Documents/disruption-py/")
import numpy as np

shot_id = 1150805012
tree_name = "efit18"

efit_cols = {
    "beta_n": r'\efit_aeqdsk:betan',
    "beta_p": r'\efit_aeqdsk:betap',
    "kappa": r'\efit_aeqdsk:eout',
    "li": r'\efit_aeqdsk:li',
    "upper_gap": r'\efit_aeqdsk:otop',
    "lower_gap": r'\efit_aeqdsk:obott',
    "q0": r'\efit_aeqdsk:q0',
    "qstar": r'\efit_aeqdsk:qstar',
    "q95": r'\efit_aeqdsk:q95',
    "v_loop_efit": r'\efit_aeqdsk:vloopt',
    "Wmhd": r'\efit_aeqdsk:wplasm',
    "ssep": r'\efit_aeqdsk:ssep',
    "n_over_ncrit": r'\efit_aeqdsk:xnnc',
    "tritop": r'\efit_aeqdsk:doutu',
    "tribot":  r'\efit_aeqdsk:doutl',
    "a_minor": r'\efit_aeqdsk:aminor',
    "rmagx":r'\efit_aeqdsk:rmagx', #TODO: change units to [m] (current [cm])
    "chisq":r'\efit_aeqdsk:chisq'
}

# Old Client
from MDSplus import Tree

efit_tree = Tree(tree=tree_name, shot=shot_id)
efit_time = efit_tree.getNode(r'\efit_aeqdsk:time').data().astype('float64', copy=False) # [s]
efit_data = dict()

#Get data from each of the columns in efit_cols one at a time
for param in efit_cols:
    try:
        efit_data[param] = efit_tree.getNode(efit_cols[param]).data().astype('float64', copy=False)
    except Exception as e:
        print(f"Failed to get efit data for param {param} with error {e}")

print(f"Retrieved {len(efit_data)} params with old client")

from MDSplus import Connection
conn = Connection('alcdata-new')
from MDSplus import Tree
import MDSplus as mds


conn.openTree("efit18", shot=shot_id)
efit_time = conn.get(r'\efit_aeqdsk:time', None).data().astype('float64', copy=False) # [s]
efit_data = dict()

#Get data from each of the columns in efit_cols one at a time
for param in efit_cols:
    try:
        efit_data[param] = conn.get(efit_cols[param]).data().astype('float64', copy=False)
    except Exception as e:
        print(f"Failed to get efit data for param {param} with error {e}")
  
print(f"Retrieved {len(efit_data)} params with thin client")

conn.openTree("pcs", shot=shot_id)
root_nid = conn.get('GetDefaultNid()')
print(root_nid)

children_nids = conn.get('getnci(getnci($, "CHILDREN_NIDS"), "NID_NUMBER")', root_nid)
print(children_nids)

desired_segment_nids, desired_segment_paths = [], []
children_paths = conn.get('getnci($, "FULLPATH")', children_nids)
for child_nid, node_path in zip(children_nids, children_paths):
    node_path = node_path.strip()
    if node_path.split(".")[-1].startswith("SEG_"):
        desired_segment_nids.append(child_nid)
        desired_segment_paths.append(node_path)

print(desired_segment_paths)

segments_is_on = conn.get(f'getnci($, "STATE")', np.array(desired_segment_nids)).data()
print(segments_is_on)
