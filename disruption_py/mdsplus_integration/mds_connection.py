from typing import List
import MDSplus as mds

class MDSConnection:
    
    def __init__(self, conn_string : str, shot_id : int, manage_tree_opening : bool = True):
        self.conn = mds.Connection(conn_string)
        self.shot_id = shot_id
        self.last_open_tree = None
        self.manage_tree_opening = manage_tree_opening
        
    def open_tree(self, tree_name : str):
        """
        Open the specified tree
        """
        if self.manage_tree_opening:
            if self.last_open_tree != tree_name:
                self.conn.openTree(tree_name, self.shot_id)
        else:
            if not self.manage_tree_opening:
                self.conn.openTree(tree_name, self.shot_id)
            
        self.last_open_tree = tree_name
    
    def close_tree(self, tree_name : str):
        """
        Close the specified tree
        """
        if self.manage_tree_opening:
            if self.last_open_tree != tree_name:
                self.conn.closeTree(tree_name, self.shot_id)
        else:
            self.conn.closeTree(tree_name, self.shot_id)
        
        self.last_open_tree = None
            
        
    def close_all_trees(self):
        """
        Close all open trees
        """
        self.last_open_tree = None
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
        if tree_name is not None:
            self.open_tree(tree_name)
        return self.conn.get(expression, arguments)