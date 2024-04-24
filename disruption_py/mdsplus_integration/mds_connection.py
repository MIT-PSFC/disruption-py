import logging
from typing import Any, Callable, Dict, List, Tuple
import MDSplus as mds
import h5pyd as h5py
import random

class ProcessMDSConnection():
    """
    Abstract class for connecting to MDSplus.
    
    Ensure that a single MDSPlus connection is used by each process for all shots retrieved by that process.
    """
    
    def __init__(self, conn_string : str):
        self.conn = None
        if conn_string != 'DoNotConnect':
            self.conn = mds.Connection(conn_string)
            
    def get_shot_connection(self, shot_id : int):
        """ Get MDSPlus Connection wrapper for individual shot. """
        return MDSConnection(self.conn, shot_id)

class HDF:
    """
    Class to fill hsds cache with answers
    """
    def __init__(self):
        self.file = None
        self.shot_id = None
        self.endpoints = [ 'http://localhost:5101', 'http://localhost:5102' ]

    def add_cache(self, tree, shot, expression, args, value):
        import h5pyd as h5py
        import random
        if self.shot_id != shot:
            self.file = h5py.File(f'/cmod/{shot}', 'a', use_cache=False, end_point=random.choice(self.endpoints))
            self.shot_id = shot
        if tree not in self.file:
            root = self.file.create_group(tree)
        else:
            root = self.file[tree]
        if args is None:
            key = f'{expression}'
        else:
            key = f'{expression}_{args}'
        if key in root:
            del root[key]
        root.create_dataset(key, data=value)
     
class MDSConnection:
    """ 
    Wrapper class for MDSPlus Connection class used for handling individual shots. 
    """
    
    logger = logging.getLogger('disruption_py')
    
    def __init__(self, conn : mds.Connection, shot_id : int):
        self.conn = conn
        self.shot_id = shot_id
        self.tree_nickname_funcs = {}
        self.tree_nicknames = {}
        self.open_trees = set()
        self.last_open_tree = None
        self.fill_hsds = False 
        self.fill_mongodb = False 
        self.use_hsds = False 
        self.use_mongo = False
        self.cache_miss_enable = False 
#        self.endpoints = [ 'http://localhost:5101', 'http://localhost:5102' ]
        self.endpoints = ['http://a00b48dc37a994baaad23bc7e01aefba-924539795.us-east-1.elb.amazonaws.com/']
        self.hdf = HDF()
 
    def open_tree(self, tree_name : str):
        """
        Open the specified _name.
        
        If the specified tree_name is nickname for a tree_name, will open the tree
        that it is a nickname for.
        """
        if tree_name not in self.tree_nicknames and tree_name in self.tree_nickname_funcs:
            self.tree_nicknames[tree_name] = self.tree_nickname_funcs[tree_name]()
        if tree_name in self.tree_nicknames:
            print("two")
            tree_name = self.tree_nicknames[tree_name]
        if not self.use_hsds and not self.use_mongo or self.cache_miss_enable:
            self.conn.openTree(tree_name, self.shot_id)
        if self.use_hsds:
            try:
                self.hdf.file = h5py.File(f'/cmod/{self.shot_id}', 'r', use_cache=False, end_point=random.choice(self.endpoints))
            except Exception as e:
                print(e)
        elif self.use_mongo:
            pass
    

        if self.conn and self.last_open_tree != tree_name:
            self.conn.openTree(tree_name, self.shot_id)
            
        self.last_open_tree = tree_name
        self.open_trees.add(tree_name)
    
    def close_tree(self, tree_name : str):
        """
        Close the specified tree_name.
        """
        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]
        
        if tree_name in self.open_trees:
            if not self.use_hsds and not self.use_mongo or self.cache_miss_enable:
                try:
                    self.conn.closeTree(tree_name, self.shot_id)
                except Exception as e:
                    self.logger.warning(f"Error closing tree {tree_name} in shot {self.shot_id}")
                    self.logger.debug(e)
                
        if self.last_open_tree == tree_name:
            self.last_open_tree = None
        self.open_trees.discard(tree_name)
        
    def close_all_trees(self):
        """
        Close all open trees
        """
        for open_tree in self.open_trees.copy():
            self.close_tree(open_tree)
        self.last_open_tree = None
        self.open_trees.clear()
        # self.conn.closeAllTrees()
        
      
    def _get(self, expression : str, arguments : Any = None, tree_name : str = None) -> Any:
        """Evaluate the specified expression.
        
        The expression is passed as string argument, but may contain optional arguments. 
        These arguments are then passed as an array of Data objects.

        Parameters
        ----------
        expression : str
            MDSplus TDI expression. Please see MDSplus documentation for more information.
        arguments : Any, optional
            Arguments for MDSplus TDI Expression. Please see MDSplus documentation for more information. Default None.
        tree_name : str, optional

        Returns
        -------
        Any
            Result of evaluating TDI expression from MDSplus.
        """
        if tree_name is not None:
            self.open_tree(tree_name)
        if not self.use_hsds and not self.use_mongo:
            ans = self.conn.get(expression, arguments).data()
        if self.use_hsds:
            try:
                if expression.startswith('_sig='):
                    expression = expression[6:]
                    self.sig_expression = expression
                elif expression.startswith('dim_of(_sig)'):
                    expression = f"dim_of({self.sig_expression}{expression[12:]}"
                ans = self.hdf.file[tree_name][f'{expression}'].value
            except:
                try:
                    if self.cache_miss_enable:
                        ans = self.conn.get(expression, arguments).data()
                    else:
                        ans = None
                except:
                    ans = None
        elif self.use_mongo:
            ans = None
        return ans
    
    # Added Methods
    def get(self, expression : str, arguments : Any = None, tree_name : str = None) -> Any:
        """Evaluate the specified expression.
        
        The expression is passed as string argument, but may contain optional arguments. 
        These arguments are then passed as an array of Data objects.

        Parameters
        ----------
        expression : str
            MDSplus TDI expression. Please see MDSplus documentation for more information.
        arguments : Any, optional
            Arguments for MDSplus TDI Expression. Please see MDSplus documentation for more information. Default None.
        tree_name : str, optional

        Returns
        -------
        Any
            Result of evaluating TDI expression from MDSplus.
        """

        ans = self._get(expression, arguments, tree_name)
        if ans is not None:
            if self.fill_hsds:
                self.hdf.add_cache(tree_name, self.shot_id, expression, arguments, ans)
            if self.fill_mongodb:
                pass        
        return ans
    
    def get_record_data(self, path : str, tree_name : str = None, dim_nums : List = None) -> Tuple:
        """Get data and dimension for record at specified path

        Parameters
        ----------
        path : str
            MDSplus path to record.
        tree_name : str, optional
            The name of the tree that must be open for retrieval.
        dim_nums : List, optional
            A list of dimensions that should have their size retrieved. Default [0].
            
        Returns
        -------
        Tuple
            Returns the node data, followed by the requested dimensions
        """
        dim_nums = dim_nums or [0]
        
        data = self._get("_sig=" + path, treename=tree_name)
        dims = [self._get(f"dim_of(_sig,{dim_num})", treename=tree_name) for dim_num in dim_nums]
        if self.fill_hsds:
            self.hdf.add_cache(tree_name, self.shot_id, path, None, data)
        if self.fill_mongodb:
            pass
        if self.fill_hsds:
            for count,d in enumerate(dims):
                self.hdf.add_cache(tree_name, self.shot_id, f'dim_of({path}, {count})', None, d)
        if self.fill_mongodb:
            pass
        return data, *dims
    
    
    def get_dims(self, path : str, tree_name : str = None, dim_nums : List = None) -> Tuple:
        """Get the size of specified dimensions for record at specified path
        Parameters
        ----------
        path : str
            MDSplus path to record.
        tree_name : str, optional
            The name of the tree that must be open for retrieval.
        dim_nums : List, optional
            A list of dimensions that should have their size retrieved. Default [0].

        Returns
        -------
        Tuple
            Returns the requested dimensions as a tuple.
        """
        dim_nums = dim_nums or [0]

        if tree_name is not None:
            self.open_tree(tree_name)
        dims = [self.conn._get(f"dim_of({path},{dim_num})").data() for dim_num in dim_nums]
        if self.fill_hsds:
            for count,d in enumerate(dims):
                self.hdf.add_cache(tree_name, self.shot_id, f'dim_of({path}, {count})', None, d)
        if self.fill_mongodb:
            pass
        return dims
    
    # MDSplusML cache population and use
    def set_aux_cache_parameters(self, fill_hsds : bool, fill_mongodb : bool, use_hsds : bool, use_mongo : bool, cache_miss_enable : bool):
        self.fill_hsds = fill_hsds
        self.fill_mongodb = fill_mongodb
        self.use_hsds = use_hsds
        self.use_mongo = use_mongo
        self.cache_miss_enable = cache_miss_enable

    # nicknames
    
    def add_tree_nickname_funcs(self, tree_nickname_funcs : Dict[str, Callable]):
        """
        Add tree nickname functions to the connection.
        
        Required because some tree nickname function require the connection to exist.
        """
        self.tree_nickname_funcs.update(tree_nickname_funcs)
    
    def get_tree_name_of_nickname(self, nickname : str):
        """
        Get the tree name that the nickname has been set to or None if the nickname was not set.
        """
        if nickname not in self.tree_nicknames and nickname in self.tree_nickname_funcs:
            self.tree_nicknames[nickname] = self.tree_nickname_funcs[nickname]()

        return self.tree_nicknames.get(nickname, None)
    
    def tree_name(self, for_name: str) -> str:
        """
        The tree name for for_name, whether it is a nickname or tree name itself
        """
        return self.get_tree_name_of_nickname(for_name) or for_name
