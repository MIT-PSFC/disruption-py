from MDSplus import Tree, TreeNode
import logging
         
class TreeManager:
    logger = logging.getLogger('disruption_py')

    def __init__(self, shot_id):
        self._shot_id = int(shot_id)
        self._open_trees = {}
        self._closed_trees = set()
        self._nicknames = {}
        self.num_times_opened = {}
    
    def try_open_and_nickname(self, tree_name: str, nickname: str) -> Tree:
        
        tree = self.open_tree(tree_name)
        
        if nickname in self._nicknames:
            self.logger.info(f"Nickname {nickname} for tree {self._nicknames[nickname]} replaced with {tree_name}")
        
        self._nicknames[nickname] = tree_name
        
        return tree
    
    def tree_from_nickname(self, nickname:str) -> Tree:
        tree_name = self.tree_name_of_nickname(nickname)
        
        if tree_name is None:
            self.logger.info(f"No tree name exists for nickname: {nickname}")
            return None
        if tree_name not in self._open_trees:
            self.logger.info(f"Tree named {tree_name} for nickname {nickname} not open")
        
        return self.open_tree(tree_name)
    
    def tree_name_of_nickname(self, nickname:str) -> str:
        return self._nicknames.get(nickname, None)
    
    def open_tree(self, tree_name: str) -> Tree:
        self.num_times_opened[tree_name] = self.num_times_opened.get(tree_name, 0) + 1
        
        if tree_name in self._open_trees:
            return self._open_trees[tree_name]
        if tree_name in self._closed_trees:
            self.logger.info(f"Tree {tree_name} reopened after close for shot {self._shot_id}")
        
        try:
            self._open_trees[tree_name] = TreeWrapper(tree_manager=self, tree=tree_name, shot=self._shot_id, mode="readonly")
        except Exception as e:
            self.logger.debug(f"Failed to open tree {tree_name} | num trees open {len(self._open_trees)}")
            raise e
                    
        return self._open_trees[tree_name]
    
    # close trees early
    def cleanup_not_needed(self):
        total_use_estimate = {'analysis': 2, 'LH': 1, 'RF': 1, 'spectroscopy': 1, 'magnetics': 5, 'pcs': 2, 'electrons': 3, 'engineering': 1, 'xtomo': 1, 'hybrid': 1}
        current_use = [(key, value) for key, value in self.num_times_opened.items()]
        for tree_name, num_times_used in current_use:
            if total_use_estimate.get(tree_name, 0) <= num_times_used and tree_name not in self._nicknames.values():
                self.close_tree(tree_name)
    
    def cleanup(self):
        for tree_name in list(self._open_trees.keys()):
            self.close_tree(tree_name)
            
    def close_tree(self, tree_name: str):
        tree_to_close = self.mark_tree_closed(tree_name)
        if tree_to_close is not None:
            tree_to_close.close()
            
    
    def mark_tree_closed(self, tree_name: str):
        if tree_name in self._open_trees:
            tree_to_close = self._open_trees.pop(tree_name)
            self._closed_trees.add(tree_name)
            return tree_to_close
        return None
            
    ### Tree Nodes
    # @staticmethod
    # def get_data_for_tree_node(tree_node: TreeNode):
        
            
class TreeWrapper(Tree):
    
    def __init__(self, tree_manager: TreeManager, *args, **kwargs):
        self._tree_manager = tree_manager
        super().__init__(*args, **kwargs)
        
        
    def close(self, *args, **kwargs):
        if self.tree is not None:
            self._tree_manager.mark_tree_closed(self.tree)
            
    