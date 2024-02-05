from typing import List
import MDSplus as mds

class ShotConnection:
    
    def __init__(self, conn : mds.Connection, shot_id : int):
        self.conn = conn
        self.shot_id = shot_id
        
    def open_tree(self, tree_name : str):
        """
        Open the specified tree
        """
        self.conn.openTree(tree_name, self.shot_id)
    
    def close_tree(self, tree_name : str):
        """
        Close the specified tree
        """
        self.conn.closeTree(tree_name, self.shot_id)
        
    def close_all_trees(self):
        """
        Close all open trees
        """
        self.conn.closeAllTrees(self.shot_id)
        
    def set_default(self, path : str):
        """
        Set the default position in the currently open tree
        """
        self.conn.setDefault(path)
        
    def get(self, expression : str, arguments = None, tree_name : str = None):
        """
        Eevaluate the specified expression
        
        The expression is passed as string argument, but may contain optional arguments. 
        These arguments are then passed as an array of Data objects.
        """
        return self.conn.get(expression, arguments)