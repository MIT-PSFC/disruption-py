#!/usr/bin/env python3

"""
Module for populating shot data by executing physics methods in DisruptionPy.
"""

import time
from collections.abc import Iterable

import numpy as np
import xarray as xr
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
    columns = retrieval_settings.run_columns
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
) -> dict | xr.DataArray | xr.Dataset:
    """
    Execute a physics method and store the result.

    Parameters
    ----------
    physics_method_params : PhysicsMethodParams
        Parameters containing MDS connection and shot information
    bound_method_metadata : BoundMethodMetadata
        The metadata for a physics method like the associated tokamak, columns, etc.

    Returns
    -------
    dict | xr.DataArray | xr.Dataset
        The result of the executed method.
    """

    method = bound_method_metadata.bound_method
    name = bound_method_metadata.name
    physics_method_params.logger.trace("Starting method: {name}", name=name)

    try:

        result = method(params=physics_method_params)

    # pylint: disable-next=broad-exception-caught
    except Exception as e:

        # log exception
        level = "ERROR"
        if isinstance(e, mdsExceptions.MDSplusERROR):
            pass
        elif isinstance(e, (mdsExceptions.MdsException, CalculationError)):
            level = "WARNING"
        physics_method_params.logger.log(level, "{name}: {exc}", name=name, exc=repr(e))
        physics_method_params.logger.opt(exception=True).debug(name)

        # mock-up data
        result = {col: [np.nan] for col in bound_method_metadata.columns}

        # reconnect if needed
        if isinstance(e, mdsExceptions.MDSplusERROR):
            physics_method_params.mds_conn.reconnect()

    return result


def populate_shot(
    retrieval_settings: RetrievalSettings,
    physics_method_params: PhysicsMethodParams,
) -> xr.Dataset:
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
    xr.Dataset
        A dataset containing the queried data.
    """

    # concatenate built-in and user-provided classes/methods
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
    datasets = []
    for bound_method_metadata in run_bound_method_metadata:

        # run method
        result = populate_method(
            physics_method_params=physics_method_params,
            bound_method_metadata=bound_method_metadata,
        )

        # convert non-dataset dict to dataset
        if not isinstance(result, (xr.DataArray, xr.Dataset)):

            times = physics_method_params.times

            # create data_vars dict
            data = {}
            for k, v in result.items():
                if len(v) == len(times):
                    # as expected
                    data[k] = v
                    continue
                if all(np.isnan(v)):
                    # pad all-nan var
                    physics_method_params.logger.debug("All-nan data: {col}", col=k)
                    data[k] = np.nan * times
                    continue
                physics_method_params.logger.warning(
                    "Data length mismatch: {col}", col=k
                )

            # create dataset
            result = physics_method_params.to_dataset(data=data)

        datasets += [result]

    # merge dataarrays/datasets into dataset
    dataset = xr.merge(datasets)

    # log statistics
    if dataset:
        tot = len(dataset.data_vars)
        nok = int(dataset.notnull().any().to_array().sum())
    else:
        nok, tot = 0, 0
    percent = (nok / tot * 100) if tot else 0
    if percent >= 75:
        level = "SUCCESS"
        quant = "all" if percent == 100 else "most"
    elif percent >= 25:
        level = "WARNING"
        quant = "many" if percent >= 50 else "some"
    else:
        level = "ERROR"
        quant = "no" if percent == 0 else "few"

    physics_method_params.logger.log(
        level,
        "{level}! {quant} parameters have data: {nok}/{tot} ({percent:.2f}%)"
        " in {elapsed}",
        level=level.capitalize(),
        quant=quant.capitalize(),
        nok=nok,
        tot=tot,
        percent=percent,
        elapsed=get_elapsed_time(time.time() - start_time),
    )

    # select columns
    if retrieval_settings.only_requested_columns and retrieval_settings.run_columns:
        dataset = dataset[retrieval_settings.run_columns]

    # sort columns
    return dataset[sorted(dataset)]
