#!/usr/bin/env python3

"""
Module for populating shot data by executing physics methods in DisruptionPy.
"""

import time
import traceback
from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd
from MDSplus import mdsExceptions

from disruption_py.config import config
from disruption_py.core.physics_method.caching import manually_cache
from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.metadata import (
    BoundMethodMetadata,
    get_method_metadata,
    is_parametered_method,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.machine.method_holders import get_method_holders
from disruption_py.settings.retrieval_settings import RetrievalSettings

REQUIRED_COLS = {"time", "shot", "commit_hash"}


def get_prefilled_shot_data(physics_method_params: PhysicsMethodParams) -> pd.DataFrame:
    """
    Retrieve pre-filled shot data for the given physics method parameters.

    Parameters
    ----------
    physics_method_params : PhysicsMethodParams
        Parameters containing MDS connection and shot information

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the pre-filled shot data.
    """
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
            "[Shot %s]: Computation on different timebase than pre-filled shot data",
            physics_method_params.shot_id,
        )
    return pre_filled_shot_data


def get_all_physics_methods(all_passed: list) -> set:
    """
    Retrieve all callable physics methods from the provided list.

    Parameters
    ----------
    all_passed : list
        A list of objects to search for callable physics methods.

    Returns
    -------
    set
        A set of callable physics methods found in the provided list.
    """
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
) -> list:
    """
    Bind metadata to the provided physics methods.

    Parameters
    ----------
    physics_methods : set
        A set of callable physics methods to bind metadata to.
    physics_method_params : PhysicsMethodParams
        The parameters to be passed to the methods.

    Returns
    -------
    list
        A list of bound method metadata instances.
    """
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
) -> list:
    """
    Filter the bound methods to determine which should be executed.

    Parameters
    ----------
    all_bound_method_metadata : list[BoundMethodMetadata]
        A list of bound method metadata instances.
    retrieval_settings : RetrievalSettings
        The settings that dictate which methods should be run based on tags, methods,
        and columns.
    physics_method_params : PhysicsMethodParams
        The parameters that will be passed to the methods that are run.

    Returns
    -------
    list
        A list of bound method metadata instances that are eligible to run.
    """
    tags = retrieval_settings.run_tags
    methods = retrieval_settings.run_methods
    columns = REQUIRED_COLS.union(retrieval_settings.run_columns)

    methods_to_run = []
    for bound_method_metadata in all_bound_method_metadata:
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
                "[Shot %s]: Skipping %s in class %s",
                physics_method_params.shot_id,
                bound_method_metadata.name,
                bound_method_metadata.bound_method,
            )
    return methods_to_run


def populate_method(
    physics_method_params: PhysicsMethodParams,
    bound_method_metadata: BoundMethodMetadata,
    start_time: float,
) -> Any:
    """
    Execute a physics method and log the results.

    Parameters
    ----------
    physics_method_params : PhysicsMethodParams
        Parameters containing MDS connection and shot information
    bound_method_metadata : BoundMethodMetadata
        The metadata for a physics method like the associated tokamak, columns, etc.
    start_time : float
        The start time for measuring execution duration.

    Returns
    -------
    Any
        The result of the executed method, or None if an error occurred.
    """
    method = bound_method_metadata.bound_method
    name = bound_method_metadata.name

    result = None
    physics_method_params.logger.info(
        "[Shot %s]: Populating %s", physics_method_params.shot_id, name
    )
    try:
        result = method(params=physics_method_params)
    except (
        mdsExceptions.TreeNNF,
        mdsExceptions.TreeNODATA,
        CalculationError,
        NotImplementedError,
        ValueError,
    ) as e:
        if isinstance(e, ValueError):
            catch_error_msgs = [
                "x and y arrays must be equal in length along interpolation axis."
            ]
            if not any(msg in str(e.args) for msg in catch_error_msgs):
                raise

        physics_method_params.logger.warning(
            "[Shot %s]: Failed to populate %s with error %s",
            physics_method_params.shot_id,
            name,
            e,
        )
        physics_method_params.logger.debug("%s", traceback.format_exc())
        result = {col: [np.nan] for col in bound_method_metadata.columns}

    physics_method_params.logger.info(
        "[Shot %s]: Completed %s, time_elapsed: %s",
        physics_method_params.shot_id,
        name,
        time.time() - start_time,
    )
    return result


def populate_shot(
    retrieval_settings: RetrievalSettings,
    physics_method_params: PhysicsMethodParams,
) -> pd.DataFrame:
    """
    Run the physics methods to populate shot data.

    This function executes the physics methods included through the
    `custom_physics_methods` property of retrieval_settings or in the built-in list
    of methods. It selects methods based on run_methods, run_tags, and run_columns
    in retrieval_settings.

    Parameters
    ----------
    retrieval_settings : RetrievalSettings
        The shot settings dictating what methods should be run.
    physics_method_params : PhysicsMethodParams
        Parameter that will be passed to methods that are run.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the queried data.
    """
    # Concatenate built-in classes containing registered methods with user-provided
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
    # Manually cache data that has already been retrieved (likely from SQL tables)
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
                        "[Shot %s]: Skipping %s already populated",
                        physics_method_params.shot_id,
                        method_metadata.name,
                    )

    start_time = time.time()
    methods_data = []
    for bound_method_metadata in run_bound_method_metadata:
        if bound_method_metadata in cached_method_metadata:
            continue
        methods_data.append(
            populate_method(
                physics_method_params=physics_method_params,
                bound_method_metadata=bound_method_metadata,
                start_time=start_time,
            )
        )

    filtered_methods = []
    for method_dict in methods_data:
        if method_dict is None:
            continue
        # Pad parameters which are only a single NaN (from our error outputs) in
        # order to create a DataFrame for easy comparison with cached data.
        for parameter in method_dict:
            if (
                np.all(np.isnan(method_dict[parameter]))
                and len(method_dict[parameter]) == 1
            ):
                method_dict[parameter] = np.full(len(pre_filled_shot_data), np.nan)
        method_df = pd.DataFrame(method_dict)
        if len(method_df) != len(pre_filled_shot_data):
            physics_method_params.logger.error(
                "[Shot %s]: Ignoring parameter %s with different length than timebase",
                physics_method_params.shot_id,
                method_dict,
            )
            # TODO: Should we drop the columns, or is it better to raise an
            # exception when the data do not match?
            continue
        filtered_methods.append(method_df)

    # TODO: This is a hack to get around the fact that some methods return
    #       multiple parameters. This should be fixed in the future.

    local_data = pd.concat([pre_filled_shot_data] + filtered_methods, axis=1)
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
