#!/usr/bin/env python3

"""
This module provides a decorator to signify methods that calculate physics quantities.
"""

import time
from functools import wraps
from typing import Callable, List, Union

from disruption_py.core.physics_method.caching import cache_method
from disruption_py.core.physics_method.metadata import MethodMetadata
from disruption_py.machine.tokamak import Tokamak


def physics_method(
    cache: bool = True,
    tokamak: Union[Tokamak, List[Tokamak]] = None,
    columns: Union[List[str], Callable] = None,
) -> Callable:
    """
    Decorator to signify a method to be run by DisruptionPy.

    The decorated method calculates disruption parameters and returns a Dataset.
    All decorated methods must take the single argument params of type
    `PhysicsMethodParams`. The decorated method will be run if designated by the
    `run_methods` or `run_columns` attributes of the `RetrievalSettings`
    class, and if included inside of the `custom_physics_methods` argument of the
    `retrieval_settings` or in the built-in method list. If run, the result of the
    decorated method will be output to the `output_setting`.

    A common pattern for parameterized methods is first retrieving data from MDSplus
    using the `mds_conn` and then using that retrieved data to compute data to return.

    Parameters
    ----------
    cache : bool, optional
        Whether to cache the result of the method, by default True.
    tokamak : Union['Tokamak', List['Tokamak']], optional
        A list of Tokamak objects that represent which tokamak this method may
        be used for, by default None, allows the method to be run for any tokamak.
    columns : Union[List[str], Callable], optional
        The columns that are in the DataFrame returned by the method. Alternatively,
        you can pass a method that returns the names of used trees at runtime. Default
        value is an empty list implying that no columns are returned by the function.

    Returns
    -------
    Callable
        A decorated method with caching and metadata attributes.
    """

    def outer_wrapper(method: Callable) -> Callable:
        if cache:
            wrapper = cache_method(method)
        else:
            wrapper = method

        # Add elapsed time log
        @wraps(wrapper)
        def timed_wrapper(params, *args, **kwargs):
            start_time = time.time()
            result = wrapper(params, *args, **kwargs)
            params.logger.verbose(
                "{t:.3f}s : {name}",
                name=method.__name__,
                t=time.time() - start_time,
            )
            return result

        method_metadata = MethodMetadata(
            name=method.__name__,
            cache=cache,
            tokamaks=tokamak,
            columns=columns,
        )

        timed_wrapper.method_metadata = method_metadata

        return timed_wrapper

    return outer_wrapper
