from disruption_py.database import create_cmod_handler
import sqlite3
import pandas as pd

data_handler = create_cmod_handler()
shot_table = data_handler.query(f"select * from shots", use_pandas=True)
data_handler.get_shot(shot_table.iloc[0]["SHOT"])


# fetch all SQL tables:
conn = data_handler.conn
cursor = conn.cursor()
cursor.execute("SELECT table_name FROM information_schema.tables")
print(cursor.fetchall())

def get_shot_data_from_sql_table(self, shot_ids=None, cols=None, sql_table=None):
    shot_ids = ','.join([str(shot_id) for shot_id in shot_ids])
    cols = ' '.join([f"{col}," for col in cols[:-1]]) + f" {cols[-1]}"
    if cols is None:
        cols = "*"
    if shot_ids is None:
        query = f"select {cols} from {sql_table}"
    else:
        query = f"select {cols} from {sql_table} where shot in ({shot_ids})"
    shot_df = pd.read_sql_query(query, self.conn)
    return shot_df


from disruption_py.database import create_d3d_handler   
data_handler = create_d3d_handler()
shot_table = data_handler.query(f"select * from shots", use_pandas=True)

data_handler.get_efit_tree(195749)
data_handler.get_efit_trees(shot_table.iloc[0:100]["SHOT"].values)


# There's no table called plasmas, so pass on this code. Let's assume 'EFIT01' is the only efit tree. 
# conn = data_handler.conn
# with conn.cursor() as curs:
#     curs.execute(f"select tree from plasmas where shot in (40607, 40608, 40609) and deleted = 0 order by idx")
#     efit_trees = curs.fetchall()

# NOTE: data_handler.get_shot() only queries disruption_warning.

test_shot = d3d_shot.D3DShot(shot_id=195749, efit_tree_name='EFIT01')
test_shot.data

# ^-- that one was empty. Trying again. 

test_shot = d3d_shot.D3DShot(shot_id=190586, efit_tree_name='EFIT01')
test_shot.data

# ^-- same 

test_shot = d3d_shot.D3DShot(shot_id=180586, efit_tree_name='EFIT01')
test_shot.data

# use the generate_dataset.py function. 
# Collate a list of sample test shots in disruption_warning

disruption_warning_shots = data_handler.query("select * from disruption_warning where shot in (1150403001)", use_pandas=True)


# possible tree names: efit01, 'd3d', efitrt1

import MDSplus

## D3D
shot_number = 156201
tree_name = "EFIT01"
tree = MDSplus.Tree(tree_name, shot_number)

## CMOD
shot_number = 1160330014  
tree_name = "EFIT01" 

# Open the MDSPlus tree
tree = MDSplus.Tree(tree_name, shot_number)

# Get the top-level node for the shot
top_node = tree.getNode('\SHOT')

# Get all the descendants (parameters) of the top-level node
descendants = top_node.getDescendants()

# Iterate over the descendants and print their names
for param in descendants:
    print(param.getNodeName())

# ^-- above does not work.

tree.getNodeWild("***") ## -- prints out a list of numbers, unevenly up to 1681!!

import matplotlib.pyplot as plt

wild_nodes = tree.getNodeWild("***")
wild_node_names = []
for node in wild_nodes:
    print(node, tree.getNode(node))
    wild_node_names.append(tree.getNode(node))

# convert wild_node_names to a numpy array and save as text
import numpy as np
np.array(wild_node_names).tofile("wild_node_names.txt", sep="\n")



# Tests: make sure that interpolation doesn't shift the dataset
# Zander is starting w generate_dataset.py
#   VSCode on CMOD machines
#   Lucas -- share my shot numbers and columns with Zander to compare with his data (from generate_dataset.py), from Cmod_shot_class

# double check that select operations are OK with SQL 



