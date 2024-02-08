from typing import Callable, Dict, List
import MDSplus as mds

class ProcessMDSConnection():
    """
    Abstract class for connecting to MDSplus.
    
    Remove need to import mdsplus in handlers.
    """
    
    def __init__(self, conn_string : str):
        self.conn = mds.Connection(conn_string)
    
    def get_shot_connection(self, shot_id : int):
        return MDSConnection(self.conn, shot_id)
    
class MDSConnection:
    
    def __init__(self, conn : mds.Connection, shot_id : int):
        self.conn = conn
        self.shot_id = shot_id
        self.tree_nickname_funcs = {}
        self.tree_nicknames = {}
        self.last_open_tree = None    
    
    def add_tree_nickname_funcs(self, tree_nickname_funcs : Dict[str, Callable]):
        """
        Add tree nickname functions to the connection.
        
        Required because some tree nickname function require the connection to exist.
        """
        self.tree_nickname_funcs.update(tree_nickname_funcs)
    
    def open_tree(self, tree_name : str):
        """
        Open the specified tree
        """
        if tree_name not in self.tree_nicknames and tree_name in self.tree_nickname_funcs:
            self.tree_nicknames[tree_name] = self.tree_nickname_funcs[tree_name]()
            
        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]
        
        if self.last_open_tree != tree_name:
            self.conn.openTree(tree_name, self.shot_id)
            
        self.last_open_tree = tree_name
    
    def close_tree(self, tree_name : str):
        """
        Close the specified tree
        """
        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]
            
        if self.last_open_tree == tree_name:
            self.last_open_tree = None
        self.conn.closeTree(tree_name, self.shot_id)
        
    def close_all_trees(self):
        """
        Close all open trees
        """
        self.last_open_tree = None
        # self.conn.closeAllTrees()
        
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