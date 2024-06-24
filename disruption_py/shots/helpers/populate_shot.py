#!/usr/bin/env python3

from collections.abc import Iterable
import time
import traceback

import numpy as np
import pandas as pd

from disruption_py.settings.shot_data_request import (
    ShotDataRequestParams,
)
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.shots.helpers.method_metadata import (
    BoundMethodMetadata,
    get_method_metadata,
    is_registered_method,
)
from disruption_py.shots.helpers.method_caching import manually_cache
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


def get_prefilled_shot_data(shot_props: ShotProps):
    pre_filled_shot_data = shot_props.pre_filled_shot_data

    # If the shot object was already passed data in the constructor, use that data. Otherwise, create an empty dataframe.
    if pre_filled_shot_data is None:
        pre_filled_shot_data = pd.DataFrame()
    if "time" not in pre_filled_shot_data:
        pre_filled_shot_data["time"] = shot_props.times
    if "shot" not in pre_filled_shot_data:
        pre_filled_shot_data["shot"] = int(shot_props.shot_id)
    pre_filled_shot_data["commit_hash"] = shot_props.metadata.get("commit_hash", None)

    # Check that pre_filled_shot_data is on the same timebase as the shot object to ensure data consistency
    if (
        len(pre_filled_shot_data["time"]) != len(shot_props.times)
        or not np.isclose(
            pre_filled_shot_data["time"], shot_props.times, atol=TIME_CONST
        ).all()
    ):
        shot_props.logger.error(
            f"[Shot {shot_props.shot_id}]: ERROR Computation on different timebase than pre-filled shot data"
        )
    return pre_filled_shot_data


def get_all_registered_methods(all_passed: list):
    registered_methods = set()
    for passed in all_passed:
        if callable(passed) and is_registered_method(passed):
            registered_methods.add(passed)

        for method_name in dir(passed):
            method = getattr(passed, method_name, None)
            if method is None or not is_registered_method(method):
                continue
            registered_methods.add(method)
    return registered_methods


def bind_method_metadata(registered_methods: set, params: ShotDataRequestParams):
    all_bound_method_metadata = []
    for method in registered_methods:
        method_metadata = get_method_metadata(method, should_throw=True)
        bound_method_metadata = BoundMethodMetadata.bind(
            method_metadata=method_metadata,
            bound_method=method,
            params=params,
        )
        all_bound_method_metadata.append(bound_method_metadata)
    return all_bound_method_metadata


def filter_methods_to_run(
    all_bound_method_metadata: list[BoundMethodMetadata],
    shot_settings: ShotSettings,
    params: ShotDataRequestParams,
):
    shot_props = params.shot_props
    tags = shot_settings.run_tags
    methods = shot_settings.run_methods
    columns = REQUIRED_COLS.union(shot_settings.run_columns)

    methods_to_run = []
    for bound_method_metadata in all_bound_method_metadata:
        if not bound_method_metadata.populate:
            continue

        # exclude if tokamak does not match
        if not (
            bound_method_metadata.tokamaks is None
            or params.tokamak is bound_method_metadata.tokamaks
            or (
                isinstance(bound_method_metadata.tokamaks, Iterable)
                and params.tokamak in bound_method_metadata.tokamaks
            )
        ):
            continue

        if tags is not None and bool(
            set(bound_method_metadata.tags).intersection(tags)
        ):
            methods_to_run.append(bound_method_metadata)
        elif methods is not None and bound_method_metadata.name in methods:
            methods_to_run.append(bound_method_metadata)
        elif columns is not None and bool(
            set(bound_method_metadata.columns).intersection(columns)
        ):
            methods_to_run.append(bound_method_metadata)
        else:
            params.logger.info(
                f"[Shot {shot_props.shot_id}]:Skipping {bound_method_metadata.name} in class {bound_method_metadata.bound_method}"
            )
    return methods_to_run


def populate_method(
    params: ShotDataRequestParams,
    bound_method_metadata: BoundMethodMetadata,
    start_time,
):

    shot_props = params.shot_props
    method = bound_method_metadata.bound_method
    name = bound_method_metadata.name

    result = None
    if bound_method_metadata.populate:
        params.logger.info(f"[Shot {shot_props.shot_id}]:Populating {name}")
        try:
            result = method(params=params)
        except Exception as e:
            params.logger.warning(
                f"[Shot {shot_props.shot_id}]:Failed to populate {name} with error {e}"
            )
            params.logger.debug(f"{traceback.format_exc()}")
    else:
        params.logger.info(f"[Shot {shot_props.shot_id}]:Caching {name}")
        try:
            method(params=params)
        except Exception as e:
            params.logger.warning(
                f"[Shot {shot_props.shot_id}]:Failed to cache {name} with error {e}"
            )
            params.logger.debug(f"{traceback.format_exc()}")

    params.logger.info(
        f"[Shot {shot_props.shot_id}]:Completed {name}, time_elapsed: {time.time() - start_time}"
    )
    return result


def populate_shot(
    shot_settings: ShotSettings,
    params: ShotDataRequestParams,
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
    # Concatanate built in clases containing registred methods, with user provided classes/methods
    all_shot_data_request = (
        built_in_method_factory(params.tokamak) + shot_settings.shot_data_requests
    )
    all_registered_methods = get_all_registered_methods(all_shot_data_request)
    all_bound_method_metadata: list[BoundMethodMetadata] = bind_method_metadata(
        all_registered_methods, params
    )
    run_bound_method_metadata: list[BoundMethodMetadata] = filter_methods_to_run(
        all_bound_method_metadata, shot_settings, params
    )

    pre_filled_shot_data = get_prefilled_shot_data(shot_props)
    # Manually cache data that has already been retrieved (likely from sql tables)
    # Methods added to pre_cached_method_names will be skipped by method optimizer
    cached_method_metadata = []
    if shot_props.pre_filled_shot_data is not None:
        for method_metadata in all_bound_method_metadata:
            cache_success = manually_cache(
                shot_props=shot_props,
                data=pre_filled_shot_data,
                method_name=method_metadata.name,
                method_columns=method_metadata.columns,
            )
            if cache_success:
                cached_method_metadata.append(method_metadata)
                if method_metadata in run_bound_method_metadata:
                    shot_props.logger.info(
                        f"[Shot {shot_props.shot_id}]:Skipping {method_metadata.name} already populated"
                    )

    start_time = time.time()
    parameters = []
    for bound_method_metadata in run_bound_method_metadata:
        if bound_method_metadata in cached_method_metadata:
            continue
        parameters.append(
            populate_method(
                params=params,
                bound_method_metadata=bound_method_metadata,
                start_time=start_time,
            )
        )

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
