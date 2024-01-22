from collections import Counter
from MDSplus import Tree, TreeNode
import logging
from typing import Callable, Dict, List, Set, Tuple
import threading

from disruption_py.utils.environment_vars import temporary_env_vars

EnvModifications = Tuple[Tuple[str, str]]

class TreeManager:
    logger = logging.getLogger('disruption_py')

    def __init__(self, shot_id):
        self._shot_id = int(shot_id)
        self._nicknames = {}
        self._nickname_environment_modifications = {}
        
        # create tree nickname when nickname first used
        self._lazy_nickname_functions = {}
        
        self.num_times_opened = {}
        self.all_opened_tree_names = set()
        
        self._open_trees = {}
        self._closed_trees = {}
        
        # count the number of getNode of each type
        self.counter = Counter()
    
    @property
    def thread_open_trees(self):
        """
        Get the dictionary of open trees mapping tree_name to MDSplus tree for the current thread.
        
        Note: mutating this dictionary affect the tracking of the open trees by the tree manager and 
        may lead to memory leaks.
        """
        current_thread = threading.current_thread()
        return self._get_open_trees(current_thread)
    
    @property
    def thread_open_tree_names(self):
        """
        Returns a view of the open tree names
        """
        return set(self.thread_open_trees.keys())
    
    def _get_open_trees(self, thread) -> Dict:
        return self._open_trees.setdefault(thread, {})
    
    def _get_closed_trees(self, thread) -> Set:
        return self._closed_trees.setdefault(thread, set())
    
    def open_tree(self, tree_name: str) -> Tree:
        """
        Open tree with name tree_name from MDSplus, returning the opened tree.
        """
        self.num_times_opened[tree_name] = self.num_times_opened.get(tree_name, 0) + 1
        
        if tree_name in self.thread_open_trees:
            return self.thread_open_trees[tree_name]
        if tree_name in self._closed_trees:
            self.logger.info(f"Tree {tree_name} reopened after close for shot {self._shot_id}")
        
        try:
            self.thread_open_trees[tree_name] = TreeWrapper(tree_manager=self, tree=tree_name, shot=self._shot_id, mode="readonly")
        except Exception as e:
            self.logger.debug(f"Failed to open tree {tree_name} | num trees open {len(self.thread_open_trees)}")
            raise e
        
        if tree_name in self._nicknames and self._nicknames[tree_name] != tree_name:
            self._nicknames.pop(tree_name)
            self.logger.error(f"Nickname {tree_name} removed because tree opened with the same name")
        
        self.all_opened_tree_names.add(tree_name)
        return self.thread_open_trees[tree_name]
    
    def cleanup_not_needed(self, can_tree_be_closed : Callable[[str], bool]):
        """Close trees that are not expected to be used again based on method_optimizer

        Parameters
        ----------
        can_tree_be_closed : Callable[[str], bool]
            Method that returns whether the tree will be needed later based on the tree name. Closed if not needed.
        """
        for thread in list(self._open_trees.keys()):
            for tree_name in list(self._get_open_trees(thread).keys()):
                if can_tree_be_closed(tree_name):
                    self.close_tree(tree_name)
    
    def cleanup(self):
        """
        Close all open trees for all threads.
        """
        self.logger.info(f"Node usage counts: {self.counter}")
        for thread in list(self._open_trees.keys()):
            for tree_name in list(self._get_open_trees(thread).keys()):
                self.close_tree(tree_name, thread)
            self._open_trees.pop(thread, None)
            
    def close_tree(self, tree_name: str, thread=None):
        """
        Close tree with tree_name for thread.
        """
        tree_to_close = self._mark_tree_closed(tree_name, thread)
        if tree_to_close is not None:
            tree_to_close.close()
            
    def _mark_tree_closed(self, tree_name: str, thread=None):
        """
        Remove tree from open trees dictionary and add to closed trees set for thread.
        """
        thread = thread or threading.current_thread()
        if tree_name in self._get_open_trees(thread):
            tree_to_close = self._get_open_trees(thread).pop(tree_name)
            self._get_closed_trees(thread).add(tree_name)
            return tree_to_close
        return None
    
    def nickname(self, nickname: str, tree_names_to_try: List[str], env_modifications_to_try: List[EnvModifications] = None):
        """Lazily create nickname for tree.
        
        Nicknames allow for the use of trees without knowing what the exact name of the tree will be for the given shot. 
        Upon using the nickname, tests provided tree names and environments until one works and sets the nickname to be for that
        tree name with that environment.
        For each environment in env_modifications_to_try tries to open every tree in tree_names_to_try.
        
        Parameters
        ----------
        nickname : str
            The nickname that you would like to set.
        tree_names_to_try : List[str]
            The list of tree names to test.
        env_modifications_to_try : List[EnvModifications]
            The list of environment modifications to test. Default is a single environment with no modifications.

        Raises
        ------
        Exception
            If none of the provided tree names can be opened for any of the provided environments raises exception.
        """
        if env_modifications_to_try is None or env_modifications_to_try == []:
            env_modifications_to_try = [()]
            
        def lazy_nickname_func():
            
            for env_modification in env_modifications_to_try:
                with temporary_env_vars(env_modification):
                    for tree_name in tree_names_to_try:
                        try:
                            return self._try_open_and_nickname(nickname=nickname, tree_name=tree_name, environment_modifications=env_modification)
                        except Exception as e:
                            self.logger.warning(
                                f"[Shot {self._shot_id}]:Failed to open tree {tree_name} for nickname {nickname}, with error {e}")
            
            raise Exception(f"Failed to find a valid tree name for nickname {nickname}")

        self._lazy_nickname_functions[nickname] = lazy_nickname_func
        return
        
            
    def _try_open_and_nickname(self, nickname: str, tree_name: str, environment_modifications : EnvModifications = None) -> Tree:
        """_try_open_and_nickname tries to open the tree with name tree_name, setting the nickname if successful.

        Parameters
        ----------
        nickname : str
            The nickname that should be set
        tree_name : str
            The tree name that disruption_py will attempt to be opened.
        environment_modifications : EnvModifications
            The environment modifications that will be used.

        Returns
        -------
        Tree
            The opened MDSplus tree.

        Raises
        ------
        Exception
            If the nickname is the name of a tree that is already open.
            The tree name fails to open.
        """
        if nickname in self.all_opened_tree_names and tree_name != nickname:
            self.logger.error(f"Cannot hide tree_name {nickname} with that nickname for {tree_name}")
            raise Exception(f"Cannot hide tree_name {nickname} with that nickname for {tree_name}")
        
        tree = self.open_tree(tree_name)
        
        if nickname in self._nicknames and self._nicknames[nickname] != tree_name:
            self.logger.info(f"Nickname {nickname} for tree {self._nicknames[nickname]} replaced with {tree_name}")
        
        self._nicknames[nickname] = tree_name
        self._nickname_environment_modifications[nickname] = environment_modifications
        return tree
    
    def is_nickname_created(self, nickname: str) -> bool:
        """
        Returns whether a nickname has been created.
        """
        return nickname in self._nicknames or nickname in self._lazy_nickname_functions
        
    def tree_from_nickname(self, nickname:str) -> Tree:
        """
        Get the MDSplus tree for the nickname.
        """
        tree_name = self.tree_name_of_nickname(nickname)
            
        if tree_name is None:
            self.logger.debug(f"No tree name exists for nickname: {nickname} (expected on first open)")
            return None
        if tree_name not in self.all_opened_tree_names:
            self.logger.debug(f"Tree named {tree_name} for nickname {nickname} not open")
        
        if self._nickname_environment_modifications.get(nickname, None) is not None:
            with temporary_env_vars(self._nickname_environment_modifications[nickname]):
                return self.open_tree(tree_name)
        else:
            return self.open_tree(tree_name)
    
    def tree_name_of_nickname(self, nickname:str) -> str:
        """
        Get the tree name that the nickname has been set to or None if the nickname was not set.
        """
        if nickname not in self._nicknames and nickname in self._lazy_nickname_functions:
            self.logger.info(f"Lazy nickname {nickname} creation attempted")
            self._lazy_nickname_functions[nickname]()
            
        return self._nicknames.get(nickname, None)
    
    def unique_name(self, for_name: str) -> str:
        """
        The tree name for for_name, whether it is a nickname or tree name itself
        """
        return self.tree_name_of_nickname(for_name) or for_name
        
            
class TreeWrapper(Tree):
    
    def __init__(self, tree_manager: TreeManager, *args, **kwargs):
        self._tree_manager = tree_manager
        super().__init__(*args, **kwargs)
        
    
    def getNode(self, name):
        self._tree_manager.counter[(self.tree, name)] += 1
        return super().getNode(name)

    def close(self, *args, **kwargs):
        if self.tree is not None:
            self._tree_manager._mark_tree_closed(self.tree)
            
    