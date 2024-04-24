
import logging
import MDSplus as mds
import os
import random

from typing import Any, Callable, Dict, List, Tuple

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
    def __init__(self, shot_id : int, filling : bool):
        import h5pyd as h5py
        import random
        import os

        self.shot_id = shot_id

        if os.environ.get('HSDS_ENDPOINTS') is not None:
            self.endpoints = os.environ['HSDS_ENDPOINTS'].split(',')
        else:
            self.endpoints = [ 'http://mfedata-archives:5101', 'http://mfedata-archives:5102' ]

        self.mode = 'a' if filling else 'r'
        self.file = h5py.File(f'/cmod/{self.shot_id}', self.mode, use_cache=(not filling), endpoint=random.choice(self.endpoints))

    def add_cache(self, tree, expression, args, value):
        import numpy as np
        if tree not in self.file:
            root = self.file.create_group(tree)
        else:
            root = self.file[tree]
        if args is None:
            key = f'{expression}'
        else:
            key = f'{expression}_{args}'

        key.replace('/', '\\/')
        key = str(key)
        if key in root:
            del root[key]
        if value.__class__ == np.str_:
            value = str(value)
        root.create_dataset(key, data=value)
     
    def get(self, tree, expression, args):
        try:
            if args is None:
                key = f'{expression}'
            else:
                key = f'{expression}_{args}'

            key.replace('/', '\\/')

            print(tree, key)
            return self.file[tree][key].value

        except:
            return None

class Mongo:
    """
    Class to handle the MongoDB connection
    """

    def __init__(self, shot_id : int):

        import pickle
        import pymongo

        self.shot_id = shot_id

        self.mongo = pymongo.MongoClient(
            'mongodb://mfews-slwalsh2',
            username=os.environ['MONGODB_USERNAME'],
            password=os.environ['MONGODB_PASSWORD'],
        )
        self.database = self.mongo['disruption-py']
        self.cache_table = self.database['cache']

        # from datetime import datetime
        # before = datetime.now()

        query = {
            "shot": self.shot_id,
        }

        multidoc = {}

        self.cache = {}
        for doc in self.cache_table.find(query):
            if doc['count'] > 1:

                if doc['tree'] not in multidoc:
                    multidoc[doc['tree']] = {}

                if doc['expr'] not in multidoc[doc['tree']]:
                    multidoc[doc['tree']][doc['expr']] = [None] * doc['count']

                multidoc[doc['tree']][doc['expr']][doc['index']] = doc['answer']
            
            else:

                answer = pickle.loads(doc['answer'])

                if doc['tree'] not in self.cache:
                    self.cache[doc['tree']] = {}

                if answer is not None:
                    self.cache[doc['tree']][doc['expr'].lower()] = answer

        for tree, exprs in multidoc.items():
            for expr, parts in exprs.items():

                if None in parts:
                    raise Exception('Freak out')
                
                answer = b''.join(parts)
                answer = pickle.loads(answer)

                self.cache[tree][expr] = answer

        # for tree, signals in self.cache.items():
        #     for signal in signals.keys():
        #         print(tree, signal)
        
        # print()
        # print()

        # count = 0
        # for _, signals in self.cache.items():
        #     count += len(signals)

        # after = datetime.now()
        # self.logger.warning(f"Fetching {count} signals for shot {self.shot_id} took {after - before}")

    def get(self, tree, expression, args):
        tree = tree.lower()
        expression = expression.lower()

        if args is not None:
            return None

        if tree in self.cache:
            if expression in self.cache[tree]:
                return self.cache[tree][expression]
        
        print('Missing', tree, expression, args)
        return None

    def add_cache(self, tree, expression, args, value):
        import math
        import pickle
        import pymongo
        from bson.binary import Binary

        # print('Caching', tree, expression, args)

        RECORD_SIZE_LIMIT = 15_728_640 # 15 MiB, the actual max is 16 MiB

        if args is not None:
            # Mongo doesn't currently support caching expressions with arguments
            return

        # TODO: Don't call this every time
        self.cache_table.create_index([ 'shot' ])
        self.cache_table.create_index([ 'tree', 'shot' ])
        self.cache_table.create_index([ 'tree', 'shot', 'expr' ])
        self.cache_table.create_index([ 'tree', 'shot', 'expr', 'index' ], unique=True)

        answer = Binary(pickle.dumps(value, protocol=2), subtype=128)
        count = math.ceil(len(answer) / RECORD_SIZE_LIMIT)

        query = {
            "tree": tree.lower(),
            "shot": self.shot_id,
            "expr": expression.lower(),
        }

        try:

            self.cache_table.delete_many(query)

            i = 0
            for offset in range(0, len(answer), RECORD_SIZE_LIMIT):
                doc = query.copy()
                doc['index'] = i
                doc['count'] = count
                doc['answer'] = answer[offset : offset + RECORD_SIZE_LIMIT]

                self.cache_table.insert_one(doc)
                i += 1

                # print(f"Cached {doc['tree']} {doc['shot']} {doc['expr']} #{doc['index']}")

        except Exception as e:
            print(e)

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
        self.use_hsds = False
        self.use_mongo = False
        self.fill_hsds = False
        self.fill_mongo = False
        self.cache_miss_enable = False

    def open_tree(self, tree_name : str):
        """
        Open the specified _name.

        If the specified tree_name is nickname for a tree_name, will open the tree
        that it is a nickname for.
        """
        if tree_name not in self.tree_nicknames and tree_name in self.tree_nickname_funcs:
            self.tree_nicknames[tree_name] = self.tree_nickname_funcs[tree_name]()

        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]

        if self.use_mdsplus:
            if self.last_open_tree != tree_name:
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
            if self.use_mdsplus:
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
        if self.conn:
            del self.conn
            self.conn = None

        # if self.use_mdsplus:
        #     self.conn.closeAllTrees()

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

        ans = None

        if tree_name is None:
            tree_name = self.last_open_tree

        if self.use_hsds:
            ans = self.hdf.get(tree_name, expression, arguments)

        elif self.use_mongo:
            if tree_name in self.tree_nicknames:
                tree_name = self.tree_nicknames[tree_name]

            ans = self.mongo.get(tree_name, expression, arguments)

        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]

        if tree_name is not None:
            self.open_tree(tree_name)

        if ans is None and self.use_mdsplus:
            ans = self.conn.get(expression, arguments).data()

        if ans is not None:
            if self.fill_hsds:
                self.hdf.add_cache(tree_name, expression, arguments, ans)

            if self.fill_mongo:
                self.mongo.add_cache(tree_name, expression, arguments, ans)

        if ans is None:
            raise mds.TreeNNF()

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

        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]

        if tree_name is not None:
            self.open_tree(tree_name)

        if tree_name is None:
            tree_name = self.last_open_tree

        dim_nums = dim_nums or [0]
        path = path.lower()

        data = None
        dims = []

        if self.use_hsds or self.use_mongo:
            try:
                data = self.get(path)
                dims = [ self.get(f'dim_of({path}, {dim_num})') for dim_num in dim_nums]
            
            except mds.TreeNNF:
                data = None
                dims = []

        if self.use_mdsplus and (data is None or len(dims) < len(dim_nums)):
            # Avoid self.get() to avoid caching _sig
            data = self.conn.get("_sig=" + path)
            dims = [self.conn.get(f"dim_of(_sig,{dim_num})") for dim_num in dim_nums]

        if data is None or len(dims) < len(dim_nums):
            raise mds.TreeNNF()

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

        if tree_name in self.tree_nicknames:
            tree_name = self.tree_nicknames[tree_name]

        path = path.lower()

        if tree_name is not None:
            self.open_tree(tree_name)

        dims = [ self.get(f'dim_of({path}, {dim_num})') for dim_num in dim_nums]
        return dims

    # MDSplusML cache population and use
    def set_aux_cache_parameters(self, fill_hsds : bool, fill_mongo : bool, use_hsds : bool, use_mongo : bool, cache_miss_enable : bool):
        self.fill_hsds = fill_hsds
        self.fill_mongo = fill_mongo
        self.use_hsds = use_hsds
        self.use_mongo = use_mongo
        self.cache_miss_enable = cache_miss_enable

        # SLW Can this be part of __init__() ?

        # In order to use fill_hsds or fill_mongo, you need to set cache_miss_enable
        if self.fill_hsds or self.fill_mongo:
            self.cache_miss_enable = True

        self.use_mdsplus = self.cache_miss_enable or (not self.use_hsds and not self.use_mongo)

        if self.use_hsds or self.fill_hsds:
            self.hdf = HDF(self.shot_id, self.fill_hsds)

        if self.use_mongo or self.fill_mongo:
            self.mongo = Mongo(self.shot_id)

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
