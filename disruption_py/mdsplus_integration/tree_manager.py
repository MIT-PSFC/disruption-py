from MDSplus import Tree
import logging
        
class TreeManager:
    logger = logging.getLogger('disruption_py')

    def __init__(self, shot_id):
        self._shot_id = int(shot_id)
        self._open_trees = {}
        self._closed_trees = {}
        self._nicknames = {}
    
    def open_tree(self, tree_name: str, nickname=None) -> Tree:
        if tree_name in self._open_trees:
            return self._open_trees[tree_name]
        if tree_name in self._closed_trees:
            self.logger.info(f"Tree {tree_name} reopened after close for shot {self._shot_id}")
        
        self._open_trees[tree_name] = TreeWrapper(tree_manager=self, tree=tree_name, shot=self._shot_id, mode="readonly")
        
        if nickname is not None:
            if nickname in self._nicknames and self._nicknames[nickname] != tree_name:
                self.logger.info(f"Nickname {nickname} for tree {self._nicknames[nickname]} replaced with {tree_name}")
            self._nicknames[nickname] = tree_name
        return self._open_trees[tree_name]
    
    def tree_from_nickname(self, nickname:str) -> Tree:
        tree_name = self.tree_name_of_nickname(nickname)
        if tree_name is None or tree_name not in self._open_trees:
            return None
        
        return self._open_trees[tree_name]
    
    def tree_name_of_nickname(self, nickname:str) -> str:
        if nickname not in self._nicknames:
            return None
        
        return self._nicknames[nickname]
    
    def close_tree(self, tree_name: str):
        if tree_name in self._open_trees:
            self._open_trees[tree_name].close()
            self.mark_tree_closed(tree_name)
    
    def mark_tree_closed(self, tree_name: str):
        if tree_name in self._open_trees:
            self._closed_trees[tree_name] = self._open_trees.pop(tree_name)
            
    def cleanup(self):
        for tree_name in list(self._open_trees.keys()):
            self._open_trees[tree_name].close()
            self.mark_tree_closed(tree_name)
            
class TreeWrapper(Tree):
    
    def __init__(self, tree_manager: TreeManager, *args, **kwargs):
        self._tree_manager = tree_manager
        super().__init__(*args, **kwargs)
        
        
    def close(self, *args, **kwargs):
        if self.tree is not None:
            self._tree_manager.mark_tree_closed(self.tree)
        return super().close(*args, **kwargs)
            
    