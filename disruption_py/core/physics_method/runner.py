#!/usr/bin/env python3

import time
import traceback
from collections.abc import Iterable

import numpy as np
import pandas as pd

from disruption_py.config import config
from disruption_py.core.physics_method.caching import manually_cache
from disruption_py.core.physics_method.metadata import (
    BoundMethodMetadata,
    get_method_metadata,
    is_parametered_method,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.method_holders import get_method_holders
from disruption_py.settings.retrieval_settings import RetrievalSettings

REQUIRED_COLS = {"time", "shot", "commit_hash"}


def get_prefilled_shot_data(physics_method_params: PhysicsMethodParams):
    pre_filled_shot_data = physics_method_params.pre_filled_shot_data

    # If the shot object was already passed data in the constructor, use that data.
    # Otherwise, create an empty dataframe.
    if pre_filled_shot_data is None:
        pre_filled_shot_data = pd.DataFrame()
    if "time" not in pre_filled_shot_data:
        pre_filled_shot_data["time"] = physics_method_params.times
    if "shot" not in pre_filled_shot_data:
        pre_filled_shot_data["shot"] = int(physics_method_params.shot_id)
    pre_filled_shot_data["commit_hash"] = physics_method_params.metadata.get(
        "commit_hash", None
    )

    # Check that pre_filled_shot_data is on the same timebase as the shot object
    # to ensure data consistency
    if (
        len(pre_filled_shot_data["time"]) != len(physics_method_params.times)
        or not np.isclose(
            pre_filled_shot_data["time"],
            physics_method_params.times,
            atol=config().TIME_CONST,
        ).all()
    ):
        physics_method_params.logger.error(
            f"[Shot {physics_method_params.shot_id}]: ERROR Computation on different"
            + " timebase than pre-filled shot data"
        )
    return pre_filled_shot_data


def get_all_physics_methods(all_passed: list):
    physics_methods = set()
    for passed in all_passed:
        if callable(passed) and is_parametered_method(passed):
            physics_methods.add(passed)

        for method_name in dir(passed):
            method = getattr(passed, method_name, None)
            if method is None or not is_parametered_method(method):
                continue
            physics_methods.add(method)
    return physics_methods


def bind_method_metadata(
    physics_methods: set,
    physics_method_params: PhysicsMethodParams,
):
    all_bound_method_metadata = []
    for method in physics_methods:
        method_metadata = get_method_metadata(method, should_throw=True)
        bound_method_metadata = BoundMethodMetadata.bind(
            method_metadata=method_metadata,
            bound_method=method,
            physics_method_params=physics_method_params,
        )
        all_bound_method_metadata.append(bound_method_metadata)
    return all_bound_method_metadata


def filter_methods_to_run(
    all_bound_method_metadata: list[BoundMethodMetadata],
    retrieval_settings: RetrievalSettings,
    physics_method_params: PhysicsMethodParams,
):
    tags = retrieval_settings.run_tags
    methods = retrieval_settings.run_methods
    columns = REQUIRED_COLS.union(retrieval_settings.run_columns)

    methods_to_run = []
    for bound_method_metadata in all_bound_method_metadata:
        if not bound_method_metadata.populate:
            continue

        # exclude if tokamak does not match
        if not (
            bound_method_metadata.tokamaks is None
            or physics_method_params.tokamak is bound_method_metadata.tokamaks
            or (
                isinstance(bound_method_metadata.tokamaks, Iterable)
                and physics_method_params.tokamak in bound_method_metadata.tokamaks
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
            physics_method_params.logger.info(
                f"[Shot {physics_method_params.shot_id}]:Skipping "
                + f"{bound_method_metadata.name} in class "
                + f"{bound_method_metadata.bound_method}"
            )
    return methods_to_run


def populate_method(
    physics_method_params: PhysicsMethodParams,
    bound_method_metadata: BoundMethodMetadata,
    start_time,
):

    method = bound_method_metadata.bound_method
    name = bound_method_metadata.name

    result = None
    if bound_method_metadata.populate:
        physics_method_params.logger.info(
            f"[Shot {physics_method_params.shot_id}]:Populating {name}"
        )
        try:
            result = method(params=physics_method_params)
        except Exception as e:
            physics_method_params.logger.warning(
                f"[Shot {physics_method_params.shot_id}]:Failed to populate {name}"
                + f"with error {e}"
            )
            physics_method_params.logger.debug(f"{traceback.format_exc()}")
    else:
        physics_method_params.logger.info(
            f"[Shot {physics_method_params.shot_id}]:Caching {name}"
        )
        try:
            method(params=physics_method_params)
        except Exception as e:
            physics_method_params.logger.warning(
                f"[Shot {physics_method_params.shot_id}]:Failed to cache {name} "
                + f"with error {e}"
            )
            physics_method_params.logger.debug(f"{traceback.format_exc()}")

    physics_method_params.logger.info(
        f"[Shot {physics_method_params.shot_id}]:Completed {name}, time_elapsed: "
        + f"{time.time() - start_time}"
    )
    return result


def populate_shot(
    retrieval_settings: RetrievalSettings,
    physics_method_params: PhysicsMethodParams,
) -> pd.DataFrame:
    """populate_shot runs the physics methods either included through the
    `custom_physics_methods` property of retrieval_settings or in the built-in list
    of methods.

    Selects methods based on run_methods, run_tags, and run_columns in
    retrieval_settings. Methods execution is reordered to minimize tree openings
    and trees opened simultaniously.

    Parameters
    ----------
    retrieval_settings : RetrievalSettings
        The shot settings dictating what methods should be run.
    physics_method_params : PhysicsMethodParams
        Parameter that will be passed to methods that are run.

    Returns
    -------
    pd.DataFrame
        A dataframe containing the queried data.
    """
    # Concatanate built in clases containing registred methods, with user provided
    # classes/methods
    all_physics_method_holders = (
        get_method_holders(physics_method_params.tokamak)
        + retrieval_settings.custom_physics_methods
    )
    all_physics_methods = get_all_physics_methods(all_physics_method_holders)
    all_bound_method_metadata: list[BoundMethodMetadata] = bind_method_metadata(
        all_physics_methods, physics_method_params
    )
    run_bound_method_metadata: list[BoundMethodMetadata] = filter_methods_to_run(
        all_bound_method_metadata, retrieval_settings, physics_method_params
    )

    pre_filled_shot_data = get_prefilled_shot_data(physics_method_params)
    # Manually cache data that has already been retrieved (likely from sql tables)
    # Methods added to pre_cached_method_names will be skipped by method optimizer
    cached_method_metadata = []
    if physics_method_params.pre_filled_shot_data is not None:
        for method_metadata in all_bound_method_metadata:
            cache_success = manually_cache(
                physics_method_params=physics_method_params,
                data=pre_filled_shot_data,
                method=method_metadata.bound_method,
                method_name=method_metadata.name,
                method_columns=method_metadata.columns,
            )
            if cache_success:
                cached_method_metadata.append(method_metadata)
                if method_metadata in run_bound_method_metadata:
                    physics_method_params.logger.info(
                        f"[Shot {physics_method_params.shot_id}]:Skipping "
                        + f"{method_metadata.name} already populated"
                    )

    start_time = time.time()
    methods = []
    for bound_method_metadata in run_bound_method_metadata:
        if bound_method_metadata in cached_method_metadata:
            continue
        methods.append(
            populate_method(
                physics_method_params=physics_method_params,
                bound_method_metadata=bound_method_metadata,
                start_time=start_time,
            )
        )

    filtered_methods = []
    for method in methods:
        if method is None:
            continue
        # Pad all parameters that are only nans to an array in order to create
        # a DataFrame for easy comparison with cached data.
        for parameter in method:
            if np.all(np.isnan(method[parameter])):
                method[parameter] = np.full(len(pre_filled_shot_data), np.nan)
        if len(pd.DataFrame(method)) != len(pre_filled_shot_data):
            physics_method_params.logger.error(
                f"[Shot {physics_method_params.shot_id}]:Ignoring parameter "
                + f"{method} with different length than timebase"
            )
            continue
        filtered_methods.append(method)

    # TODO: This is a hack to get around the fact that some methods return
    #       multiple parameters. This should be fixed in the future.

    local_data = pd.concat(
        [pre_filled_shot_data] + [pd.DataFrame(d) for d in filtered_methods], axis=1
    )
    local_data = local_data.loc[:, ~local_data.columns.duplicated()]
    if retrieval_settings.only_requested_columns:
        include_columns = list(
            REQUIRED_COLS.union(
                set(retrieval_settings.run_columns).intersection(
                    set(local_data.columns)
                )
            )
        )
        local_data = local_data[include_columns]
    return local_data
