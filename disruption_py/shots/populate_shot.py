import traceback
import time
import pandas as pd
import numpy as np
import concurrent
from dataclasses import fields
from concurrent.futures import ThreadPoolExecutor

from disruption_py.shots.shot import Shot
from disruption_py.utils.method_optimizer import MethodOptimizer, CachedMethod
from disruption_py.utils.method_caching import CachedMethodParams, manually_cache
from disruption_py.settings.shot_run_settings import ShotRunSettings, ShotDataRequest, ShotDataRequestParams
from disruption_py.utils.constants import TIME_CONST
    
def populate_method(params: ShotDataRequestParams, cached_method : CachedMethod, method_optimizer : MethodOptimizer, start_time):
        
        shot = params.shot
        method = cached_method.method
        
        def run_method():
            if cached_method.built_in_to_shot:
                return method()
            else:
                return method(params=params)
        
        result = None
        if callable(method) and hasattr(method, 'cached_method_params'):
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
    for field in  fields(cached_method_params):
        field_name = field.name
        field_value = getattr(cached_method_params, field_name)
        if callable(field_value):
            new_cached_method_params_dict[field_name] = field_value(object_to_search, params)
        else:
            new_cached_method_params_dict[field_name] = field_value
    return CachedMethodParams(**new_cached_method_params_dict)

def get_cache_methods_from_object(object_to_search, shot_run_settings: ShotRunSettings, params: ShotDataRequestParams):
    shot = params.shot
    tags = shot_run_settings.run_tags
    methods = shot_run_settings.run_methods
    
    built_in_to_shot = isinstance(object_to_search, Shot)
    is_shot_data_request = isinstance(object_to_search, ShotDataRequest)
        
    if built_in_to_shot:
        methods_to_search = dir(object_to_search)
    elif is_shot_data_request:
        methods_to_search = object_to_search.get_request_methods_for_tokemak(params.tokemak)
    
    methods_to_evaluate : list[CachedMethod] = []
    all_cached_methods : list[CachedMethod] = []
    for method_name in methods_to_search:
        attribute_to_check = getattr(object_to_search, method_name)
        cached_method_params: CachedMethodParams = getattr(attribute_to_check, 'cached_method_params', None)
        
        if not callable(attribute_to_check) or cached_method_params is None:
            continue
        
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
            params.logger.info(
                    f"[Shot {shot.get_shot_id()}]:Skipping {method_name} in class {object_to_search.__class__.__name__}")
    return methods_to_evaluate, all_cached_methods
            
def populate_shot(shot_run_settings: ShotRunSettings, params: ShotDataRequestParams):
    """
    Internal method to populate the disruption parameters of a shot object. 
    This method is called by the constructor and should not be called directly. It loops through all methods of the Shot class and calls the ones that have a `populate` attribute set to True and satisfy the tags and methods arguments.
    """

    shot = params.shot
    populated_data = params.existing_data
    
    # If the shot object was already passed data in the constructor, use that data. Otherwise, create an empty dataframe.
    if populated_data is None:
        populated_data = pd.DataFrame() 
    if 'time' not in populated_data:
        populated_data['time'] = shot.get_times()
    if 'shot' not in populated_data:
        populated_data['shot'] = shot.get_shot_id()

    # Loop through each attribute and find methods that should populate the shot object.
    methods_to_evaluate : list[CachedMethod] = []
    all_cached_methods : list[CachedMethod] = []
    
    methods_to_evaluate, all_cached_methods = get_cache_methods_from_object(shot, shot_run_settings, params)
    
    # Add the methods gound from the passed ShotDataRequest objects
    for shot_data_request in shot_run_settings.shot_data_requests:
        req_methods_to_evaluate, req_all_cached_methods = get_cache_methods_from_object(shot_data_request, shot_run_settings, params)
        methods_to_evaluate.extend(req_methods_to_evaluate)
        all_cached_methods.extend(req_all_cached_methods)
        
    # Check that existing data is on the same timebase as the shot object to ensure data consistency
    if not np.isclose(shot.data['time'], shot.get_times(), atol=TIME_CONST).all():
        params.logger.error(f"[Shot {shot.get_shot_id()}]: ERROR Computation on different timebase than used existing data")
        
    # Manually cache data that has already been retrieved (likely from sql tables)
    # Methods added to pre_cached_method_names will be skipped by method optimizer
    pre_cached_method_names = []
    if shot.initialized_with_data:
        for cached_method in all_cached_methods:
            cache_success = manually_cache(
                shot=shot, 
                data=shot.data, 
                method_name=cached_method.name, 
                method_columns=cached_method.computed_cached_method_params.columns
            )
            if cache_success:
                pre_cached_method_names.append(cached_method.name)
                if cached_method in methods_to_evaluate:
                    shot.logger.info(
                        f"[Shot {shot.get_shot_id()}]:Skipping {cached_method.name} already populated")

    method_optimizer : MethodOptimizer = MethodOptimizer(shot.get_tree_manager(), methods_to_evaluate, all_cached_methods, pre_cached_method_names)

    if shot.multithreading:
        futures = set()
        future_method_names = {}
        def future_for_next(next_method):
            new_future = executor.submit(populate_method, params, next_method, method_optimizer, start_time)
            futures.add(new_future)
            future_method_names[new_future] = next_method.name
        
        start_time = time.time()
        available_methods_runner = method_optimizer.get_async_available_methods_runner(future_for_next)
        with ThreadPoolExecutor(max_workers=3) as executor: 
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
                            f"[Shot {shot.get_shot_id()}]:Failed to populate {method_name} with future error {e}")
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
    local_data = pd.concat(parameters + [shot.data], axis=1)
    local_data = local_data.loc[:, ~local_data.columns.duplicated()]
    shot.data = local_data
    return local_data
