from dataclasses import fields
import traceback
import time
import pandas as pd
import numpy as np
import concurrent
from concurrent.futures import ThreadPoolExecutor
from disruption_py.shots.helpers.cached_method_params import CachedMethodParams, ParameterCachedMethodParams, get_cached_method_params, is_cached_method
from disruption_py.shots.parameter_functions import DEFAULT_SHOT_DATA_REQUESTS

from disruption_py.shots.shot_props import ShotProps
from disruption_py.shots.helpers.method_optimizer import MethodOptimizer, CachedMethod
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.method_caching import manually_cache
from disruption_py.settings.shot_data_request import ShotDataRequest, ShotDataRequestParams
from disruption_py.utils.constants import MAX_THREADS_PER_SHOT, TIME_CONST

REQUIRED_COLS = {'time', 'time_until_disrupt', 'shot', 'commit_hash'}
 
def populate_method(params: ShotDataRequestParams, cached_method : CachedMethod, start_time):
        
        shot_props = params.shot_props
        method = cached_method.method
        
        result = None
        if callable(method) and is_cached_method(method):
            params.logger.info(
                f"[Shot {shot_props.shot_id}]:Populating {cached_method.name}")
            # self._tree_manager.cleanup_not_needed()
            try:
                result = method(params=params)
            except Exception as e:
                params.logger.warning(
                    f"[Shot {shot_props.shot_id}]:Failed to populate {cached_method.name} with error {e}")
                params.logger.debug(f"{traceback.format_exc()}")
        elif callable(method) and hasattr(method, 'cached'):
            params.logger.info(
                f"[Shot {shot_props.shot_id}]:Caching {cached_method.name}")
            try:
                method(params=params)
            except Exception as e:
                params.logger.warning(
                    f"[Shot {shot_props.shot_id}]:Failed to cache {cached_method.name} with error {e}")
                params.logger.debug(f"{traceback.format_exc()}")
        else:
            params.logger.warning(
                f"[Shot {shot_props.shot_id}]:Method {cached_method.name} is not callable or does not have a `populate` attribute set to True")
            return None
        
        params.logger.info(f"[Shot {shot_props.shot_id}]:Completed {cached_method.name}, time_elapsed: {time.time() - start_time}")
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
    return cached_method_params.__class__(**new_cached_method_params_dict)

def get_cached_methods_from_object(object_to_search : ShotDataRequest, shot_settings: ShotSettings, params: ShotDataRequestParams):
    shot_props = params.shot_props
    tags = shot_settings.run_tags
    methods = shot_settings.run_methods
    columns = shot_settings.run_columns
            
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
            computed_cached_method_params=computed_cached_method_params
        )
        all_cached_methods.append(cached_method)
        
        if isinstance(computed_cached_method_params, ParameterCachedMethodParams):
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
                    f"[Shot {shot_props.shot_id}]:Skipping {method_name} in class {object_to_search.__class__.__name__}")
    return methods_to_evaluate, all_cached_methods
            
def populate_shot(shot_settings: ShotSettings, params: ShotDataRequestParams):
    """
    Internal method to populate the disruption parameters of a shot object. 
    This method is called by the constructor and should not be called directly. It loops through all methods of the Shot class and calls the ones that have a `populate` attribute set to True and satisfy the tags and methods arguments.
    """

    shot_props : ShotProps = params.shot_props
    populated_data = shot_props.populated_existing_data
    
    # If the shot object was already passed data in the constructor, use that data. Otherwise, create an empty dataframe.
    if populated_data is None:
        populated_data = pd.DataFrame() 
    if 'time' not in populated_data:
        populated_data['time'] = shot_props.times
    if 'shot' not in populated_data:
        populated_data['shot'] = int(shot_props.shot_id)
    populated_data['commit_hash'] = shot_props.metadata.get("commit_hash", None)

    # Loop through each attribute and find methods that should populate the shot object.
    methods_to_evaluate : list[CachedMethod] = []
    all_cached_methods : list[CachedMethod] = []
        
    # Add the methods gound from the passed ShotDataRequest objects
    all_shot_data_request = DEFAULT_SHOT_DATA_REQUESTS + shot_settings.shot_data_requests
    for shot_data_request in all_shot_data_request:
        req_methods_to_evaluate, req_all_cached_methods = get_cached_methods_from_object(shot_data_request, shot_settings, params)
        methods_to_evaluate.extend(req_methods_to_evaluate)
        all_cached_methods.extend(req_all_cached_methods)
        
    # Check that existing data is on the same timebase as the shot object to ensure data consistency
    if len(populated_data['time']) != len(shot_props.times) or not np.isclose(populated_data['time'], shot_props.times, atol=TIME_CONST).all():
        params.logger.error(f"[Shot {shot_props.shot_id}]: ERROR Computation on different timebase than used existing data")
        
    # Manually cache data that has already been retrieved (likely from sql tables)
    # Methods added to pre_cached_method_names will be skipped by method optimizer
    pre_cached_method_names = []
    if shot_props.populated_existing_data is not None:
        for cached_method in all_cached_methods:
            cache_success = manually_cache(
                shot_props=shot_props, 
                data=populated_data, 
                method_name=cached_method.name, 
                method_columns=cached_method.computed_cached_method_params.columns
            )
            if cache_success:
                pre_cached_method_names.append(cached_method.name)
                if cached_method in methods_to_evaluate:
                    shot_props.logger.info(
                        f"[Shot {shot_props.shot_id}]:Skipping {cached_method.name} already populated")

    method_optimizer : MethodOptimizer = MethodOptimizer(shot_props.tree_manager, methods_to_evaluate, all_cached_methods, pre_cached_method_names)

    if shot_props.num_threads_per_shot > 1:
        futures = set()
        future_method_names = {}
        def future_for_next(next_method):
            new_future = executor.submit(populate_method, params, next_method, start_time)
            futures.add(new_future)
            future_method_names[new_future] = next_method.name
        
        start_time = time.time()
        available_methods_runner = method_optimizer.get_async_available_methods_runner(future_for_next)
        with ThreadPoolExecutor(max_workers=min(shot_props.num_threads_per_shot, MAX_THREADS_PER_SHOT)) as executor: 
            available_methods_runner()
            while futures:
                done, futures = concurrent.futures.wait(futures, return_when='FIRST_COMPLETED')
                for future in done:
                    try:
                        parameter_df = future.result()
                        parameters.append(parameter_df)
                        method_optimizer.method_complete(future_method_names[future])
                        shot_props.tree_manager.cleanup_not_needed(method_optimizer.can_tree_be_closed)
                    except Exception as e:
                        params.logger.warning(
                            f"[Shot {shot_props.shot_id}]:Failed to populate {future_method_names[future]} with future error {e}")
                        params.logger.debug(
                            f"[Shot {shot_props.shot_id}: {traceback.format_exc()}")
                available_methods_runner()
                
    else:
        parameters = []
        start_time = time.time()
        method_optimizer.run_methods_sync(
            lambda next_method: parameters.append(populate_method(params, next_method, start_time))
        )
    
    filtered_parameters = []
    for parameter in parameters:
        if parameter is None:
            continue
        if len(parameter) != len(populated_data):
            params.logger.warning(
                f"[Shot {shot_props.shot_id}]:Ignoring parameters {parameter.columns} with different length than timebase")
            continue
        filtered_parameters.append(parameter)

    # TODO: This is a hack to get around the fact that some methods return
    #       multiple parameters. This should be fixed in the future.
    local_data = pd.concat(filtered_parameters + [populated_data], axis=1)
    local_data = local_data.loc[:, ~local_data.columns.duplicated()]
    if shot_settings.only_requested_columns:
        include_columns = list(REQUIRED_COLS.union(set(shot_settings.run_columns).intersection(set(local_data.columns))))
        local_data = local_data[include_columns]
    return local_data
