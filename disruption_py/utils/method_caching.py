from typing import Callable, Set, List, Dict
from dataclasses import dataclass, field
from disruption_py.mdsplus_integration.tree_manager import TreeManager
import threading

@dataclass
class CachedMethod:
    name: str
    method: Callable

def parameter_cached_method(tags=["all"], columns=[], used_trees=None, contained_cached_methods=None, **kwargs):
    """
    Decorates a function as a parameter method. Parameter methods are functions that 
    calculate disruption parameters from the data in the shot.  They are called by the Shot object when
    it is instantiated.
    """
    # TODO: Figure out how to hash _times so that we can use the cache for different timebases
    def tag_wrapper(func):
        wrapper = cached_method(used_trees=used_trees, contained_cached_methods=contained_cached_methods, **kwargs)(func)
        wrapper.populate = True
        wrapper.columns = columns
        wrapper.tags = tags
        return wrapper
    return tag_wrapper


def cached_method(used_trees=None, contained_cached_methods=None, cache_between_threads=True):
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
            if cache_between_threads:
                cache_key = func.__name__ + str(len(self._times))
            else:
                current_thread_id = threading.get_ident()
                cache_key = func.__name__ + str(len(self._times)) + str(current_thread_id)
            if cache_key in self._cached_result:
                return self._cached_result[cache_key]
            else:
                result = func(self, *args, **kwargs)
                self._cached_result[cache_key] = result
                return result

        wrapper.cached = True
        wrapper.cache_between_threads = cache_between_threads
        wrapper.used_trees = used_trees
        wrapper.contained_cached_methods = contained_cached_methods
        return wrapper
    return tag_wrapper

@dataclass
class _MethodDependecy:
    dependencies: Set[str] = field(default_factory=set)
    dependent_on: Set[str] = field(default_factory=set)
    
class MethodOptimizer:
    
    def __init__(self, tree_manager : TreeManager, parameter_methods: List[CachedMethod], all_cached_methods: List[CachedMethod]):
        self._tree_manager = tree_manager
        self._all_cached_methods : Dict[str, CachedMethod] = {method.name: method for method in all_cached_methods}
        graph: dict[str, _MethodDependecy] = {}
        visited = set(self._all_cached_methods.keys())
        cached_methods_stack = parameter_methods.copy()
        while len(cached_methods_stack) != 0:
            cached_method : CachedMethod = cached_methods_stack.pop()
            graph.setdefault(cached_method.name, _MethodDependecy())
            # contained cached methods is a list of names of used cached methods
            for contained_cached_method_name in get_or_default(cached_method.method.contained_cached_methods, []):
                if contained_cached_method_name not in self._all_cached_methods:
                    continue
                if contained_cached_method_name not in visited:
                    visited.add(contained_cached_method_name)
                    cached_methods_stack.append(self._all_cached_methods[contained_cached_method_name])
                graph.setdefault(contained_cached_method_name, _MethodDependecy()).dependencies.add(cached_method.name)
                graph[cached_method.name].dependent_on.add(contained_cached_method_name)
        self._method_dependencies = graph
        
        self._tree_remaining_count = self.calculate_tree_remaining_count()
        
        self._methods_in_progress = set()
    
    def next_method(self) -> (CachedMethod, int):
        if len(self._method_dependencies) == 0:
            return None, len(self._method_dependencies)
        # get methods that are not dependent on another cached method
        # guranteed to never be empty unless there is a cycle
        methods_to_consider = [cached_method_name for cached_method_name in self._method_dependencies 
                                if len(self._method_dependencies[cached_method_name].dependent_on) == 0]
        if len(methods_to_consider) == 0:
            return None, len(self._method_dependencies)

        method_tree_counts = []
        open_tree_names = self._tree_manager.get_open_tree_names()
        for cached_method_name in methods_to_consider:
            used_trees = self.get_used_trees_by_unique_name(cached_method_name)
            remaining_open_tree_count = sum(self._tree_remaining_count[tree_name] for tree_name in used_trees)
            close_estimate = len([tree_name for tree_name in used_trees 
                                            if self._tree_remaining_count[tree_name] == 0 and tree_name in open_tree_names])
            open_estiamte = len([tree_name for tree_name in used_trees 
                                            if self._tree_remaining_count[tree_name] != 0 and tree_name not in open_tree_names])
            method_tree_counts.append((open_estiamte - close_estimate, remaining_open_tree_count, cached_method_name))
        next_method_name = min(method_tree_counts)[2]
        
        # update self._method_dependencies for method removal
        self._method_dependencies.pop(next_method_name)

        next_method = self._all_cached_methods[next_method_name]
        
        self._methods_in_progress.add(next_method_name)
        return next_method, len(self._method_dependencies)
    
    def run_methods_sync(self, method_executor:Callable[[CachedMethod], None]):
        while (True):
            next_method, _ = self.next_method()
            if (next_method == None):
                break
            method_executor(next_method)
            self.method_complete(next_method.name)
            self._tree_manager.cleanup_not_needed(self.can_tree_be_closed)
                
    def get_async_available_methods_runner(self, method_executor:Callable[[CachedMethod], None]):
        def inner():
            to_execute = []
            while (True):
                next_method, _ = self.next_method()
                if (next_method == None): 
                    break
                to_execute.append(next_method)
            for method in to_execute:
                method_executor(method)
        return inner
    
    def method_complete(self, method_name: str):
        for method_to_check in self._method_dependencies:
            # remove method that has compeleted from all methods dependent_on
            self._method_dependencies[method_to_check].dependent_on.discard(method_name)
        # update self._tree_remaining_count
        # for tree_name in self.get_used_trees_by_unique_name(method_name):
        #     self._tree_remaining_count[tree_name] -= 1
        #     if self._tree_remaining_count[tree_name] < 0:
        #         self._tree_remaining_count.pop(tree_name)
        self._methods_in_progress.remove(method_name)
        self._tree_remaining_count = self.calculate_tree_remaining_count()
    
    def can_tree_be_closed(self, tree_name: str):
        return tree_name not in self._tree_remaining_count

    def calculate_tree_remaining_count(self):
        tree_remaining_count = {}
        for cached_method_name in set(self._method_dependencies.keys()) | getattr(self, "_methods_in_progress", set()):
            for tree_name in self.get_used_trees_by_unique_name(cached_method_name):
                # count initialized to 0, as represents remaining count after chosen method completed
                cur_count = tree_remaining_count.get(tree_name, -1)
                tree_remaining_count[tree_name] = cur_count + 1
        return tree_remaining_count
    
    def get_used_trees_by_unique_name(self, method_name):
        ''' Returns value if not None, otherwise returns default '''
        used_trees = get_or_default(self._all_cached_methods[method_name].method.used_trees, [])
        return [self._tree_manager.unique_name(used_tree) for used_tree in used_trees]
    

def get_or_default(value, default):
    ''' Returns value if not None, otherwise returns default '''
    return default if value is None else value