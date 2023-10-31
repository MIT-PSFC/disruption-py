from MDSplus import Tree, TreeNode
import logging
from typing import Callable
import threading

class TreeManager:
    logger = logging.getLogger('disruption_py')

    def __init__(self, shot_id):
        self._shot_id = int(shot_id)
        self._nicknames = {}
        self.num_times_opened = {}
        self.all_opened_tree_names = set()
        
        self._open_trees = {}
        self._closed_trees = {}
    
    @property
    def thread_open_trees(self):
        current_thread = threading.current_thread()
        return self.get_open_trees(current_thread)
    
    def get_open_trees(self, thread):
        return self._open_trees.setdefault(thread, {})
    
    def get_closed_trees(self, thread):
        return self._closed_trees.setdefault(thread, set())
    
    def try_open_and_nickname(self, tree_name: str, nickname: str) -> Tree:
        
        if nickname in self.all_opened_tree_names and tree_name != nickname:
            self.logger.error(f"Cannot hide tree_name {nickname} with that nickname for {tree_name}")
            raise f"Cannot hide tree_name {nickname} with that nickname for {tree_name}"
        
        tree = self.open_tree(tree_name)
        
        if nickname in self._nicknames and self._nicknames[nickname] != tree_name:
            self.logger.info(f"Nickname {nickname} for tree {self._nicknames[nickname]} replaced with {tree_name}")
        
        self._nicknames[nickname] = tree_name
        
        return tree
    
    def tree_from_nickname(self, nickname:str) -> Tree:
        tree_name = self.tree_name_of_nickname(nickname)
        
        if tree_name is None:
            self.logger.debug(f"No tree name exists for nickname: {nickname} (expected on first open)")
            return None
        if tree_name not in self.all_opened_tree_names:
            self.logger.debug(f"Tree named {tree_name} for nickname {nickname} not open")
        
        return self.open_tree(tree_name)
    
    def tree_name_of_nickname(self, nickname:str) -> str:
        return self._nicknames.get(nickname, None)
    
    def unique_name(self, for_name: str) -> str:
        return self.tree_name_of_nickname(for_name) or for_name
    
    def open_tree(self, tree_name: str) -> Tree:
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

    def get_open_tree_names(self):
        '''
        Get a set of all open_tree names including the nicknames of open tree names
        '''
        open_tree_names = set(self.thread_open_trees.keys())
        for nickname, tree_name in self._nicknames.items():
            if tree_name in open_tree_names:
                open_tree_names.add(nickname)
        return open_tree_names
    
    def cleanup_not_needed(self, can_tree_be_closed : Callable):
        '''
        Close trees that are not expected to be used again based on method_optimizer
        '''
        for thread in list(self._open_trees.keys()):
            for tree_name in list(self.get_open_trees(thread).keys()):
                if can_tree_be_closed(tree_name):
                    self.close_tree(tree_name)
    
    def cleanup(self):
        for thread in list(self._open_trees.keys()):
            for tree_name in list(self.get_open_trees(thread).keys()):
                self.close_tree(tree_name, thread)
            self._open_trees.pop(thread)
            
    def close_tree(self, tree_name: str, thread=None):
        tree_to_close = self.mark_tree_closed(tree_name, thread)
        if tree_to_close is not None:
            tree_to_close.close()
            
    def mark_tree_closed(self, tree_name: str, thread=None):
        thread = thread or threading.current_thread()
        if tree_name in self.get_open_trees(thread):
            tree_to_close = self.get_open_trees(thread).pop(tree_name)
            self.get_closed_trees(thread).add(tree_name)
            return tree_to_close
        return None
        
            
class TreeWrapper(Tree):
    
    def __init__(self, tree_manager: TreeManager, *args, **kwargs):
        self._tree_manager = tree_manager
        super().__init__(*args, **kwargs)
        
        
    def close(self, *args, **kwargs):
        if self.tree is not None:
            self._tree_manager.mark_tree_closed(self.tree)
            
    