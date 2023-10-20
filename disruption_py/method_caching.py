
from dataclasses import dataclass

@dataclass
class CachedMethod:
    name: str
    method: function

def parameter_cached_method(tags=["all"], used_trees=None, contained_cached_methods=None):
    """
    Decorates a function as a parameter method. Parameter methods are functions that 
    calculate disruption parameters from the data in the shot.  They are called by the Shot object when
    it is instantiated.
    """
    # TODO: Figure out how to hash _times so that we can use the cache for different timebases
    def tag_wrapper(func):
        wrapper = cached_method(used_trees=used_trees, contained_cached_methods=contained_cached_methods)(func)

        wrapper.populate = True
        wrapper.tags = tags
        return wrapper
    return tag_wrapper


def cached_method(used_trees=None, contained_cached_methods=None, wait_for_callers=True):
    """
    Decorates a function as a cached method and instantiates its cache. Cached methods are functions that 
    run expensive operations on data in the shot and may be reused. The cache is used to store the results 
    of the parameter method so that it is only calculated once per shot for a given timebase. 
    
    wait_for_callers: only allow for the opened resources to be released after the calling methods have completed
    """
    # TODO: Figure out how to hash _times so that we can use the cache for different timebases
    def tag_wrapper(func):
        def wrapper(self, *args, **kwargs):
            # Create the cache if it doesn't exist
            if not hasattr(self, '_cached_result'):
                self._cached_result = {}
            cache_key = func.__name__ + str(len(self._times))
            if cache_key in self._cached_result:
                return self._cached_result[cache_key]
            else:
                result = func(self, *args, **kwargs)
                self._cached_result[cache_key] = result
                return result

        wrapper.cached = True
        wrapper.used_trees = used_trees
        wrapper.contained_cached_methods = contained_cached_methods
        return wrapper
    return tag_wrapper

@dataclass
class _MethodDependecy:
    dependencies: set = set()
    dependent_on: set = set()
    
class MethodOptimizer:
    
    def __init__(self, parameter_methods: list[CachedMethod], all_cached_methods: list[CachedMethod]):
        self._all_cached_methods : dict[str, function] = {method.name: method.method for method in all_cached_methods}
        graph: dict[str, _MethodDependecy] = {}
        visited = set(self._all_cached_methods.keys())
        cached_methods_stack = parameter_methods.copy()
        while len(cached_methods_stack) != 0:
            cached_method : CachedMethod = cached_methods_stack.pop()
            graph.setdefault(cached_method.name, _MethodDependecy())
            # contained cached methods is a list of names of used cached methods
            for contained_cached_method_name in cached_method.contained_cached_methods:
                if contained_cached_method_name not in self._all_cached_methods:
                    continue
                if contained_cached_method_name not in visited:
                    visited.add(contained_cached_method_name)
                    cached_methods_stack.append(self._all_cached_methods[contained_cached_method_name])
                graph.setdefault(contained_cached_method_name, _MethodDependecy()).dependencies.add(cached_method.name)
                graph[cached_method.name].dependent_on.add(contained_cached_method_name)
        self._method_dependencies = graph
        
        # remaining uses of tree count by tree name
        tree_remaining_count = {}
        for cached_method_name in self._method_dependencies:
            cached_method : CachedMethod = self._all_cached_methods[cached_method_name]
            for tree_name in cached_method.method.used_trees:
                # count initialized to 0, as represents remaining count after chosen method completed
                cur_count = tree_remaining_count.get(tree_name, -1)
                tree_remaining_count.set(tree_name, cur_count + 1)
        self._tree_remaining_count = tree_remaining_count
        
        self._methods_in_progress = set()
    
    def next_method(self, open_tree_names: set[str]) -> (CachedMethod, int):
        if len(self._method_dependencies) == 0:
            return None, len(self._method_dependencies)
        # get methods that are not dependent on another cached method
        # guranteed to never be empty unless there is a cycle
        methods_to_consider = [cached_method_name for cached_method_name in self._method_dependencies 
                                if len(self._method_dependencies[cached_method_name].dependent_on) == 0]
        if len(methods_to_consider):
            return None, len(self._method_dependencies)

        method_tree_counts = []
        for cached_method_name in methods_to_consider:
            used_trees = self._all_cached_methods[cached_method_name].used_trees
            remaining_open_tree_count = sum(self._tree_remaining_count[tree_name] for tree_name in used_trees)
            close_estimate = len([tree_name for tree_name in used_trees 
                                            if self._tree_remaining_count[tree_name] == 0 and tree_name in open_tree_names])
            open_estiamte = len([tree_name for tree_name in used_trees 
                                            if self._tree_remaining_count[tree_name] != 0 and tree_name not in open_tree_names])
            method_tree_counts.append(open_estiamte - close_estimate, remaining_open_tree_count, cached_method_name)
        next_method_name = min(method_tree_counts)[2]
        
        # update self._method_dependencies for method removal
        self._method_dependencies.pop(next_method_name)

        next_method = self._all_cached_methods[next_method_name]
        # update self._tree_remaining_count
        for tree_name in next_method.used_trees:
            self._tree_remaining_count[tree_name] -= 1
            if self._tree_remaining_count[tree_name] < 0:
                self._tree_remaining_count.pop(tree_name)
        
        self._methods_in_progress.add(next_method_name)
        return next_method, len(self._method_dependencies)
    
    def method_complete(self, method_name: str):
        for method_dependency_name in self._method_dependencies[method_name].dependencies:
            self._method_dependencies[method_dependency_name].dependent_on.remove(method_name)
        self._methods_in_progress.remove(method_name)
    
    def can_tree_be_closed(self, tree_name: str):
        return tree_name not in self._tree_remaining_count
