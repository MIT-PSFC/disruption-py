import time
import traceback
from dataclasses import fields
from typing import List

import numpy as np
import pandas as pd

from disruption_py.settings.shot_data_request import (
    ShotDataRequest,
    ShotDataRequestParams,
)
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.cached_method_props import (
    CachedMethodParams,
    CachedMethodProps,
    ParameterCachedMethodParams,
    get_cached_method_params,
    is_cached_method,
)
from disruption_py.shots.helpers.method_caching import manually_cache
from disruption_py.shots.helpers.method_optimizer import MethodOptimizer
from disruption_py.shots.shot_props import ShotProps
from disruption_py.utils.constants import TIME_CONST
from disruption_py.utils.mappings.tokamak import Tokamak

REQUIRED_COLS = {"time", "shot", "commit_hash"}


def built_in_method_factory(tokamak: Tokamak):
    if tokamak is Tokamak.D3D:
        from disruption_py.shots.parameter_methods.d3d.built_in import (
            D3D_DEFAULT_SHOT_DATA_REQUESTS,
        )

        return D3D_DEFAULT_SHOT_DATA_REQUESTS
    elif tokamak is Tokamak.CMOD:
        from disruption_py.shots.parameter_methods.cmod.built_in import (
            CMOD_DEFAULT_SHOT_DATA_REQUESTS,
        )

        return CMOD_DEFAULT_SHOT_DATA_REQUESTS
    else:
        raise ValueError(f"Invalid tokamak for built-ins {tokamak}")


def populate_method(
    params: ShotDataRequestParams, cached_method_props: CachedMethodProps, start_time
):

    shot_props = params.shot_props
    method = cached_method_props.method

    result = None
    if callable(method) and is_cached_method(method):
        params.logger.info(
            f"[Shot {shot_props.shot_id}]:Populating {cached_method_props.name}"
        )
        try:
            result = method(params=params)
        except Exception as e:
            params.logger.warning(
                f"[Shot {shot_props.shot_id}]:Failed to populate {cached_method_props.name} with error {e}"
            )
            params.logger.debug(f"{traceback.format_exc()}")
    elif callable(method) and hasattr(method, "cached"):
        params.logger.info(
            f"[Shot {shot_props.shot_id}]:Caching {cached_method_props.name}"
        )
        try:
            method(params=params)
        except Exception as e:
            params.logger.warning(
                f"[Shot {shot_props.shot_id}]:Failed to cache {cached_method_props.name} with error {e}"
            )
            params.logger.debug(f"{traceback.format_exc()}")
    else:
        params.logger.warning(
            f"[Shot {shot_props.shot_id}]:Method {cached_method_props.name} is not callable or does not have a `populate` attribute set to True"
        )
        return None

    params.logger.info(
        f"[Shot {shot_props.shot_id}]:Completed {cached_method_props.name}, time_elapsed: {time.time() - start_time}"
    )
    return result


def compute_cached_method_params(
    cached_method_params: CachedMethodParams,
    object_to_search: ShotDataRequest,
    params: ShotDataRequestParams,
) -> CachedMethodParams:
    """Evaluate arguments to decorators to usable values.

    Some parameters provided to the cached_method and parameter_cached_method decorators can take method that are evaluated
    at runtime. `compute_cached_method_params` evaluates all of these methods and returns a new instance of `CachedMethodParams`
    without functiohn parameters.

    Parameters
    ----------
    cached_method_params : CachedMethodParams
        The parameters that we check for function values, and evaluate them if they exist.
    object_to_search : ShotDataRequest
        The ShotDataRequest object that contained the decorated method.
    params : ShotDataRequestParams
        Params passed to the fumnction for evaluation. These are the same parameters passed when retrieving data using decorated method.

    Returns
    -------
    CachedMethodParams
        The new instance of cached method params with all of its values having been evaluated.
    """
    new_cached_method_params_dict = {}
    for field in fields(cached_method_params):
        field_name = field.name
        field_value = getattr(cached_method_params, field_name)
        if callable(field_value):
            new_cached_method_params_dict[field_name] = field_value(
                object_to_search, params
            )
        else:
            new_cached_method_params_dict[field_name] = field_value
    return cached_method_params.__class__(**new_cached_method_params_dict)


def _get_cached_methods_from_object(
    object_to_search: ShotDataRequest,
    shot_settings: ShotSettings,
    params: ShotDataRequestParams,
):
    """Get methods decorated with cached_method or parameter_cached_method decorator for object, that should be run.

    Parameters
    ----------
    object_to_search : ShotDataRequest
        The object that hsa its properties searched for decorated methods.
    shot_settings : ShotSettings
        The shot settings dictating what methods should be run.
    params : ShotDataRequestParams
        Parameter that will be passed to methods that are run.

    Returns
    -------
    Tuple[List[CachedMethod], List[CachedMethod]]
        A list of paaneter methods that require evaluation, and a list of all decorated methods
    """
    shot_props = params.shot_props
    tags = shot_settings.run_tags
    methods = shot_settings.run_methods
    columns = REQUIRED_COLS.union(shot_settings.run_columns)

    methods_to_search = object_to_search.get_request_methods_for_tokamak(params.tokamak)

    cached_methods_to_evaluate_props: List[CachedMethodProps] = []
    all_cached_methods_props: List[CachedMethodProps] = []
    for method_name in methods_to_search:
        attribute_to_check = getattr(object_to_search, method_name)

        if not callable(attribute_to_check) or not is_cached_method(attribute_to_check):
            continue

        cached_method_params: CachedMethodParams = get_cached_method_params(
            attribute_to_check, should_throw=True
        )
        computed_cached_method_params = compute_cached_method_params(
            cached_method_params, object_to_search, params
        )
        cached_method_props = CachedMethodProps(
            name=method_name,
            method=attribute_to_check,
            computed_cached_method_params=computed_cached_method_params,
        )
        all_cached_methods_props.append(cached_method_props)

        if isinstance(computed_cached_method_params, ParameterCachedMethodParams):
            # If method does not have tag included and name included then skip
            if tags is not None and bool(
                set(computed_cached_method_params.tags).intersection(tags)
            ):
                cached_methods_to_evaluate_props.append(cached_method_props)
                continue
            if methods is not None and method_name in methods:
                cached_methods_to_evaluate_props.append(cached_method_props)
                continue
            if columns is not None and bool(
                set(computed_cached_method_params.columns).intersection(columns)
            ):
                cached_methods_to_evaluate_props.append(cached_method_props)
                continue
            params.logger.info(
                f"[Shot {shot_props.shot_id}]:Skipping {method_name} in class {object_to_search.__class__.__name__}"
            )
    return cached_methods_to_evaluate_props, all_cached_methods_props


def populate_shot(
    shot_settings: ShotSettings, params: ShotDataRequestParams
) -> pd.DataFrame:
    """populate_shot runs the parameter methods in the shot_data_requests property of shot_settings.

    Selects methods based on run_mdethods, run_tags, and run_columns in shot_settings.
    Methods execution is reordered to minimize tree openings and trees opened simultaniously.

    Parameters
    ----------
    shot_settings : ShotSettings
        The shot settings dictating what methods should be run.
    params : ShotDataRequestParams
        Parameter that will be passed to methods that are run.

    Returns
    -------
    pd.DataFrame
        A dataframe containing the querried data.
    """
    shot_props: ShotProps = params.shot_props
    pre_filled_shot_data = shot_props.pre_filled_shot_data

    # If the shot object was already passed data in the constructor, use that data. Otherwise, create an empty dataframe.
    if pre_filled_shot_data is None:
        pre_filled_shot_data = pd.DataFrame()
    if "time" not in pre_filled_shot_data:
        pre_filled_shot_data["time"] = shot_props.times
    if "shot" not in pre_filled_shot_data:
        pre_filled_shot_data["shot"] = int(shot_props.shot_id)
    pre_filled_shot_data["commit_hash"] = shot_props.metadata.get("commit_hash", None)

    # Loop through each attribute and find methods that should populate the shot object.
    cached_methods_to_evaluate_props: list[CachedMethodProps] = []
    all_cached_methods_props: list[CachedMethodProps] = []

    # Add the methods gound from the passed ShotDataRequest objects
    all_shot_data_request = (
        built_in_method_factory(params.tokamak) + shot_settings.shot_data_requests
    )
    for shot_data_request in all_shot_data_request:
        req_methods_to_evaluate, req_all_cached_methods = (
            _get_cached_methods_from_object(shot_data_request, shot_settings, params)
        )
        cached_methods_to_evaluate_props.extend(req_methods_to_evaluate)
        all_cached_methods_props.extend(req_all_cached_methods)

    # Check that pre_filled_shot_data is on the same timebase as the shot object to ensure data consistency
    if (
        len(pre_filled_shot_data["time"]) != len(shot_props.times)
        or not np.isclose(
            pre_filled_shot_data["time"], shot_props.times, atol=TIME_CONST
        ).all()
    ):
        params.logger.error(
            f"[Shot {shot_props.shot_id}]: ERROR Computation on different timebase than pre-filled shot data"
        )

    # Manually cache data that has already been retrieved (likely from sql tables)
    # Methods added to pre_cached_method_names will be skipped by method optimizer
    pre_cached_method_names = []
    if shot_props.pre_filled_shot_data is not None:
        for cached_method_props in all_cached_methods_props:
            cache_success = manually_cache(
                shot_props=shot_props,
                data=pre_filled_shot_data,
                method_name=cached_method_props.name,
                method_columns=cached_method_props.get_param_value("columns", None),
            )
            if cache_success:
                pre_cached_method_names.append(cached_method_props.name)
                if cached_method_props in cached_methods_to_evaluate_props:
                    shot_props.logger.info(
                        f"[Shot {shot_props.shot_id}]:Skipping {cached_method_props.name} already populated"
                    )

    method_optimizer: MethodOptimizer = MethodOptimizer(
        mds_conn=shot_props.mds_conn,
        parameter_cached_method_props=cached_methods_to_evaluate_props,
        all_cached_method_props=all_cached_methods_props,
        pre_cached_method_names=pre_cached_method_names,
    )

    parameters = []
    start_time = time.time()

    def next_method_runner(next_method: CachedMethodProps):
        if isinstance(
            next_method.computed_cached_method_params, ParameterCachedMethodParams
        ):
            parameters.append(populate_method(params, next_method, start_time))
        else:
            # if methods are cached_methods (meaning not parameter methods) we don't return their data
            populate_method(params, next_method, start_time)

    method_optimizer.run_methods_sync(next_method_runner)

    filtered_parameters = []
    for parameter in parameters:
        if parameter is None:
            continue
        if len(parameter) != len(pre_filled_shot_data):
            params.logger.error(
                f"[Shot {shot_props.shot_id}]:Ignoring parameter {parameter} with different length than timebase"
            )
            continue
        filtered_parameters.append(parameter)

    # TODO: This is a hack to get around the fact that some methods return
    #       multiple parameters. This should be fixed in the future.
    local_data = pd.concat([pre_filled_shot_data] + filtered_parameters, axis=1)
    local_data = local_data.loc[:, ~local_data.columns.duplicated()]
    if shot_settings.only_requested_columns:
        include_columns = list(
            REQUIRED_COLS.union(
                set(shot_settings.run_columns).intersection(set(local_data.columns))
            )
        )
        local_data = local_data[include_columns]
    return local_data
