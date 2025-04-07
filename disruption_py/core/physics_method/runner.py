#!/usr/bin/env python3

"""
Module for populating shot data by executing physics methods in DisruptionPy.
"""

import time
from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd
from MDSplus import mdsExceptions

from disruption_py.core.physics_method.errors import CalculationError
from disruption_py.core.physics_method.metadata import (
    BoundMethodMetadata,
    get_method_metadata,
    is_physics_method,
)
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.misc import get_elapsed_time
from disruption_py.machine.method_holders import get_method_holders
from disruption_py.settings.retrieval_settings import RetrievalSettings

REQUIRED_COLS = {"shot", "time"}


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
        if callable(passed) and is_physics_method(passed):
            physics_methods.add(passed)

        for method_name in dir(passed):
            method = getattr(passed, method_name, None)
            if method is None or not is_physics_method(method):
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
        The settings that dictate which methods should be run based on methods
        and columns.
    physics_method_params : PhysicsMethodParams
        The parameters that will be passed to the methods that are run.

    Returns
    -------
    list
        A list of bound method metadata instances that are eligible to run.
    """
    methods = retrieval_settings.run_methods
    if retrieval_settings.run_columns is not None:
        columns = REQUIRED_COLS.union(retrieval_settings.run_columns)
    else:
        columns = None
    only_excluded_methods_specified = all("~" in method for method in methods or [])
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

        both_none = methods is None and columns is None
        method_specified = methods is not None and bound_method_metadata.name in methods
        column_specified = columns is not None and bool(
            set(bound_method_metadata.columns).intersection(columns)
        )
        is_not_excluded = (
            only_excluded_methods_specified
            and not columns
            and methods
            and f"~{bound_method_metadata.name}" not in methods
        )
        should_run = (
            both_none or method_specified or column_specified or is_not_excluded
        )

        # reasons that methods should be excluded from should run
        should_not_run = (
            methods is not None and f"~{bound_method_metadata.name}" in methods
        )

        if should_run and not should_not_run:
            methods_to_run.append(bound_method_metadata)
        else:
            physics_method_params.logger.debug(
                "Skipping method: {name}", name=bound_method_metadata.name
            )
    return methods_to_run


def populate_method(
    physics_method_params: PhysicsMethodParams,
    bound_method_metadata: BoundMethodMetadata,
) -> Any:
    """
    Execute a physics method and log the results.

    Parameters
    ----------
    physics_method_params : PhysicsMethodParams
        Parameters containing MDS connection and shot information
    bound_method_metadata : BoundMethodMetadata
        The metadata for a physics method like the associated tokamak, columns, etc.

    Returns
    -------
    Any
        The result of the executed method, or None if an error occurred.
    """
    method = bound_method_metadata.bound_method
    name = bound_method_metadata.name

    physics_method_params.logger.trace("Starting method: {name}", name=name)

    try:
        result = method(params=physics_method_params)
    # pylint: disable-next=broad-exception-caught
    except Exception as e:
        level = "ERROR"
        if isinstance(e, (mdsExceptions.MdsException, CalculationError)):
            if not isinstance(e, mdsExceptions.MDSplusERROR):
                level = "WARNING"
        physics_method_params.logger.log(level, "{name}: {exc}", name=name, exc=repr(e))
        physics_method_params.logger.opt(exception=True).debug(name)
        result = {col: [np.nan] for col in bound_method_metadata.columns}

    return result


def populate_shot(
    retrieval_settings: RetrievalSettings,
    physics_method_params: PhysicsMethodParams,
) -> pd.DataFrame:
    """
    Run the physics methods to populate shot data.

    This function executes the physics methods included through the
    `custom_physics_methods` property of retrieval_settings or in the built-in list
    of methods. It selects methods based on run_methods and run_columns
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

    # run methods and collect data
    start_time = time.time()
    methods_data = [
        populate_method(
            physics_method_params=physics_method_params,
            bound_method_metadata=bound_method_metadata,
        )
        for bound_method_metadata in run_bound_method_metadata
    ]

    # create DataFrames of proper shape
    num_parameters = 0
    num_valid = 0
    filtered_methods = []
    for method_dict in methods_data:
        if method_dict is None:
            continue
        # Pad parameters which are only a single NaN (from our error outputs) in
        # order to create a DataFrame for easy comparison with cached data.
        for parameter in method_dict:
            num_parameters += 1
            if (
                np.all(np.isnan(method_dict[parameter]))
                and len(method_dict[parameter]) == 1
            ):
                method_dict[parameter] = physics_method_params.times * np.nan
            else:
                num_valid += 1
        method_df = pd.DataFrame(method_dict)
        if len(method_df) != len(physics_method_params.times):
            physics_method_params.logger.error(
                "Ignoring parameters {parameter} with different length than timebase",
                parameter=list(method_dict.keys()),
            )
            continue
        filtered_methods.append(method_df)

    # log statistics
    percent_valid = (num_valid / num_parameters * 100) if num_parameters else 0
    if percent_valid >= 75:
        level = "SUCCESS"
        quant = "all" if percent_valid == 100 else "most"
    elif percent_valid >= 25:
        level = "WARNING"
        quant = "many" if percent_valid >= 50 else "some"
    else:
        level = "ERROR"
        quant = "no" if percent_valid == 0 else "few"

    physics_method_params.logger.log(
        level,
        "{level}! {quant} parameters have data: {num_valid}/{total} ({percent_valid:.2f}%)"
        " in {elapsed}",
        level=level.capitalize(),
        quant=quant.capitalize(),
        num_valid=num_valid,
        total=num_parameters,
        percent_valid=percent_valid,
        elapsed=get_elapsed_time(time.time() - start_time),
    )

    # concatenate partial DataFrames
    coords = pd.DataFrame(
        {
            "shot": [physics_method_params.shot_id] * len(physics_method_params.times),
            "time": physics_method_params.times,
        }
    )
    local_data = pd.concat([coords] + filtered_methods, axis=1)
    local_data = local_data.loc[:, ~local_data.columns.duplicated()]

    # include requested columns
    include_cols = set(local_data.columns).difference(REQUIRED_COLS)
    if (
        retrieval_settings.only_requested_columns
        and retrieval_settings.run_columns is not None
    ):
        include_cols = set(retrieval_settings.run_columns).intersection(include_cols)

    # sort columns
    local_data = local_data[list(REQUIRED_COLS) + sorted(list(include_cols))]
    return local_data
