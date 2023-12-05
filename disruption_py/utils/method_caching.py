from disruption_py.utils.mappings.tokemak import Tokemak
import threading
import pandas as pd
from typing import List, Callable, Union
from dataclasses import dataclass

@dataclass(frozen=True)
class CachedMethodParams:
    cache_between_threads: bool
    used_trees : Union[List[str], Callable]
    contained_cached_methods : Union[List[str], Callable]
    tokemaks : List[Tokemak]    
    # Only for parameter cached methods
    populate : bool = False
    columns : Union[List[str], Callable] = None
    tags : List[str] = None
    
    def make_parameter_method(self, columns, tags):
        return CachedMethodParams(
            cache_between_threads=self.cache_between_threads,
            used_trees=self.used_trees, 
            contained_cached_methods=self.contained_cached_methods, 
            tokemaks=self.tokemaks,
            populate=True,
            columns=columns,
            tags=tags,
        )
        
def get_method_cache_key(method_name, times):
    current_thread_id = threading.get_ident()
    return method_name + str(len(times)) + str(current_thread_id)

def parameter_cached_method(tags=["all"], columns=[], used_trees=None, contained_cached_methods=None, tokemaks=None, **kwargs):
    """
    Decorates a function as a parameter method. Parameter methods are functions that 
    calculate disruption parameters from the data in the shot.  They are called by the Shot object when
    it is instantiated.
    """
    # TODO: Figure out how to hash _times so that we can use the cache for different timebases
    def tag_wrapper(func):
        wrapper = cached_method(used_trees=used_trees, contained_cached_methods=contained_cached_methods, **kwargs)(func)
        wrapper.cached_method_params = wrapper.cached_method_params.make_parameter_method(columns, tags)
        
        return wrapper
    return tag_wrapper


def cached_method(used_trees=None, contained_cached_methods=None, cache_between_threads=True, tokemaks=None):
    """
    Decorates a function as a cached method and instantiates its cache. Cached methods are functions that 
    run expensive operations on data in the shot and may be reused. The cache is used to store the results 
    of the parameter method so that it is only calculated once per shot for a given timebase. 
    
    wait_for_callers: only allow for the opened resources to be released after the calling methods have completed
    built_in: the function is a built-in method of the Shot class and should not have shot_data_request_parameters passed as an argument
    """
    # TODO: Figure out how to hash _times so that we can use the cache for different timebases
    def tag_wrapper(func):
                
        def wrapper(self, *args, **kwargs):
            
            if "params" in kwargs:
                shot_data_request_params = kwargs["params"]
                shot = shot_data_request_params.shot
            else:
                shot_data_request_params = None
                shot = self
            
            # Create the cache if it doesn't exist
            if not hasattr(shot, '_cached_result'):
                shot._cached_result = {}
            cache_key = get_method_cache_key(func.__name__, shot._times)
            if cache_key in shot._cached_result:
                return shot._cached_result[cache_key]
            else:
                result = func(self, *args, **kwargs)
                shot._cached_result[cache_key] = result
                return result
        
        cached_method_params = CachedMethodParams(
            cache_between_threads=cache_between_threads, 
            used_trees=used_trees, 
            contained_cached_methods=contained_cached_methods,
            tokemaks=tokemaks, 
        )
        wrapper.cached_method_params = cached_method_params
        return wrapper
    return tag_wrapper


def manually_cache(shot, data : pd.DataFrame, method_name, method_columns : List[str]) -> bool:
    if method_columns is None:
        return False
    if not hasattr(shot, '_cached_result'):
        shot._cached_result = {}
    missing_columns = set(col for col in method_columns if col not in data.columns)
    if len(missing_columns) == 0:
        cache_key = get_method_cache_key(method_name, data['time'])
        shot._cached_result[cache_key] = data[method_columns]
        shot.logger.debug(
                f"[Shot {shot.get_shot_id()}]:Manually caching {method_name}")
        return True
    else:
        shot.logger.debug(
                f"[Shot {shot.get_shot_id()}]:Can not cache {method_name} missing columns {missing_columns}")
        return False