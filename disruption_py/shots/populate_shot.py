from dataclasses import fields
import traceback
import time
import pandas as pd
import numpy as np
import concurrent
from concurrent.futures import ThreadPoolExecutor

from disruption_py.shots.shot import Shot
from disruption_py.utils.method_optimizer import MethodOptimizer, CachedMethod
from disruption_py.utils.method_caching import CachedMethodParams, manually_cache, is_cached_method, get_cached_method_params
from disruption_py.settings import ShotSettings, ShotDataRequest, ShotDataRequestParams
from disruption_py.utils.constants import MAX_THREADS_PER_SHOT, TIME_CONST

REQUIRED_COLS = {'time', 'time_until_disrupt', 'shot', 'commit_hash'}
 
def populate_method(params: ShotDataRequestParams, cached_method : CachedMethod, start_time):
        
        shot = params.shot
        method = cached_method.method
        
        def run_method():
            if cached_method.built_in_to_shot:
                return method()
            else:
                return method(params=params)
        
        result = None
        if callable(method) and is_cached_method(method):
            params.logger.info(
                f"[Shot {shot.get_shot_id()}]:Populating {cached_method.name}")
            # self._tree_manager.cleanup_not_needed()
            try:
                result = run_method()
            except Exception as e:
                params.logger.warning(
                    f"[Shot {shot.get_shot_id()}]:Failed to populate {cached_method.name} with error {e}")
                params.logger.debug(f"{traceback.format_exc()}")
        elif callable(method) and hasattr(method, 'cached'):
            params.logger.info(
                f"[Shot {shot.get_shot_id()}]:Caching {cached_method.name}")
            try:
                run_method()
            except Exception as e:
                params.logger.warning(
                    f"[Shot {shot.get_shot_id()}]:Failed to cache {cached_method.name} with error {e}")
                params.logger.debug(f"{traceback.format_exc()}")
        else:
            params.logger.warning(
                f"[Shot {shot.get_shot_id()}]:Method {cached_method.name} is not callable or does not have a `populate` attribute set to True")
            return None
        
        params.logger.info(f"[Shot {shot.get_shot_id()}]:Completed {cached_method.name}, time_elapsed: {time.time() - start_time}")
        return result


def compute_cached_method_params(cached_method_params : CachedMethodParams, object_to_search, params: ShotDataRequestParams) -> CachedMethodParams:
    new_cached_method_params_dict = {}
    for field in fields(cached_method_params):
        field_name = field.name
        field_value = getattr(cached_method_params, field_name)
        if callable(field_value):
            new_cached_method_params_dict[field_name] = field_value(object_to_search, params)
        else:
            new_cached_method_params_dict[field_name] = field_value
    return CachedMethodParams(**new_cached_method_params_dict)

def get_cached_methods_from_object(object_to_search, shot_settings: ShotSettings, params: ShotDataRequestParams):
    shot = params.shot
    tags = shot_settings.run_tags
    methods = shot_settings.run_methods
    columns = shot_settings.run_columns
    
    built_in_to_shot = isinstance(object_to_search, Shot)
    is_shot_data_request = isinstance(object_to_search, ShotDataRequest)
        
    if built_in_to_shot:
        methods_to_search = dir(object_to_search)
    elif is_shot_data_request:
        methods_to_search = object_to_search.get_request_methods_for_tokamak(params.tokamak)
    
    methods_to_evaluate : list[CachedMethod] = []
    all_cached_methods : list[CachedMethod] = []
    for method_name in methods_to_search:
        attribute_to_check = getattr(object_to_search, method_name)
        
        if not callable(attribute_to_check) or not is_cached_method(attribute_to_check):
            continue
        
        cached_method_params: CachedMethodParams = get_cached_method_params(attribute_to_check, should_throw=True)
        computed_cached_method_params =  compute_cached_method_params(cached_method_params, object_to_search, params)
        cached_method = CachedMethod(
            name=method_name, 
            method=attribute_to_check, 
            built_in_to_shot=built_in_to_shot,
            computed_cached_method_params=computed_cached_method_params
        )
        all_cached_methods.append(cached_method)
        
        if computed_cached_method_params.populate is True:
            # If method does not have tag included and name included then skip
            if tags is not None and bool(set(computed_cached_method_params.tags).intersection(tags)):
                methods_to_evaluate.append(cached_method)
                continue
            if methods is not None and method_name in methods:
                methods_to_evaluate.append(cached_method)
                continue
            if columns is not None and bool(set(computed_cached_method_params.columns).intersection(columns)):
                methods_to_evaluate.append(cached_method)
                continue
            params.logger.info(
                    f"[Shot {shot.get_shot_id()}]:Skipping {method_name} in class {object_to_search.__class__.__name__}")
    return methods_to_evaluate, all_cached_methods
            
def populate_shot(shot_settings: ShotSettings, params: ShotDataRequestParams):
    """
    Internal method to populate the disruption parameters of a shot object. 
    This method is called by the constructor and should not be called directly. It loops through all methods of the Shot class and calls the ones that have a `populate` attribute set to True and satisfy the tags and methods arguments.
    """

    shot : Shot = params.shot
    populated_data = shot.data
    
    # If the shot object was already passed data in the constructor, use that data. Otherwise, create an empty dataframe.
    if populated_data is None:
        populated_data = pd.DataFrame() 
    if 'time' not in populated_data:
        populated_data['time'] = shot.get_times()
    if 'shot' not in populated_data:
        populated_data['shot'] = shot.get_shot_id()
    populated_data['commit_hash'] = shot.get_commit_hash()

    # Loop through each attribute and find methods that should populate the shot object.
    methods_to_evaluate : list[CachedMethod] = []
    all_cached_methods : list[CachedMethod] = []
    
    methods_to_evaluate, all_cached_methods = get_cached_methods_from_object(shot, shot_settings, params)
    
    # Add the methods gound from the passed ShotDataRequest objects
    for shot_data_request in shot_settings.shot_data_requests:
        req_methods_to_evaluate, req_all_cached_methods = get_cached_methods_from_object(shot_data_request, shot_settings, params)
        methods_to_evaluate.extend(req_methods_to_evaluate)
        all_cached_methods.extend(req_all_cached_methods)
        
    # Check that existing data is on the same timebase as the shot object to ensure data consistency
    if len(populated_data['time']) != len(shot.get_times()) or not np.isclose(populated_data['time'], shot.get_times(), atol=TIME_CONST).all():
        params.logger.error(f"[Shot {shot.get_shot_id()}]: ERROR Computation on different timebase than used existing data")
        
    # Manually cache data that has already been retrieved (likely from sql tables)
    # Methods added to pre_cached_method_names will be skipped by method optimizer
    pre_cached_method_names = []
    if shot.initialized_with_data:
        for cached_method in all_cached_methods:
            cache_success = manually_cache(
                shot=shot, 
                data=populated_data, 
                method_name=cached_method.name, 
                method_columns=cached_method.computed_cached_method_params.columns
            )
            if cache_success:
                pre_cached_method_names.append(cached_method.name)
                if cached_method in methods_to_evaluate:
                    shot.logger.info(
                        f"[Shot {shot.get_shot_id()}]:Skipping {cached_method.name} already populated")

    method_optimizer : MethodOptimizer = MethodOptimizer(shot.get_tree_manager(), methods_to_evaluate, all_cached_methods, pre_cached_method_names)

    if shot.num_threads_per_shot > 1:
        futures = set()
        future_method_names = {}
        def future_for_next(next_method):
            new_future = executor.submit(populate_method, params, next_method, start_time)
            futures.add(new_future)
            future_method_names[new_future] = next_method.name
        
        start_time = time.time()
        available_methods_runner = method_optimizer.get_async_available_methods_runner(future_for_next)
        with ThreadPoolExecutor(max_workers=min(shot.num_threads_per_shot, MAX_THREADS_PER_SHOT)) as executor: 
            available_methods_runner()
            while futures:
                done, futures = concurrent.futures.wait(futures, return_when='FIRST_COMPLETED')
                for future in done:
                    try:
                        parameter_df = future.result()
                        parameters.append(parameter_df)
                        method_optimizer.method_complete(future_method_names[future])
                        shot.get_tree_manager().cleanup_not_needed(method_optimizer.can_tree_be_closed)
                    except Exception as e:
                        params.logger.warning(
                            f"[Shot {shot.get_shot_id()}]:Failed to populate {future_method_names[future]} with future error {e}")
                        params.logger.debug(
                            f"[Shot {shot.get_shot_id()}: {traceback.format_exc()}")
                available_methods_runner()
                
    else:
        parameters = []
        start_time = time.time()
        method_optimizer.run_methods_sync(
            lambda next_method: parameters.append(populate_method(params, next_method, method_optimizer, start_time))
        )
                
    parameters = [
        parameter for parameter in parameters if parameter is not None]
    # TODO: This is a hack to get around the fact that some methods return
    #       multiple parameters. This should be fixed in the future.
    local_data = pd.concat(parameters + [populated_data], axis=1)
    local_data = local_data.loc[:, ~local_data.columns.duplicated()]
    if shot_settings.only_requested_columns:
        include_columns = REQUIRED_COLS.union(set(shot_settings.run_columns).intersection(set(local_data.columns)))
        local_data = local_data[include_columns]
    shot.data = local_data
    return local_data
