#!/usr/bin/env python3

"""
Module for utility functions related to class instantiation, data manipulation, and version control.
"""

import os
import subprocess
import sys
import time
import warnings
from functools import lru_cache
from pathlib import Path
from tempfile import mkdtemp
from typing import Dict, List, Tuple, Type

import numpy as np
from loguru import logger


def without_duplicates(lst: List):
    """
    Get list without duplicates while maintaining order.

    Parameters
    ----------
    lst : List
        List to get without duplicates.

    Returns
    -------
    List
        The list lst with duplicates removed.
    """
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]


def safe_cast(array: np.ndarray, dtype, copy=False) -> np.ndarray:
    """
    Safely cast a NumPy array to a specified dtype, suppressing warnings.

    Parameters
    ----------
    array : np.ndarray
        The NumPy array to cast.
    dtype : type
        The target data type to cast the array to.
    copy : bool, optional
        Whether to create a copy of the array (default is False).

    Returns
    -------
    np.ndarray
        The casted NumPy array.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        return array.astype(dtype, copy=copy)


@lru_cache
def get_commit_hash() -> str:
    """
    Retrieve the current Git commit hash.

    Returns
    -------
    str
        The commit hash, if available.
    """
    try:
        commit_hash = (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .decode("ascii")
            .strip()
        )
    except subprocess.CalledProcessError:
        commit_hash = ""
    return commit_hash


@lru_cache
def get_temporary_folder() -> str:
    """
    Create and return a temporary folder.
    The result is cached to return the same path for different invocations.

    Returns
    -------
    str
        Resulting temporary folder.
    """

    # create temporary top folder
    top = os.path.join(
        "/tmp",
        os.getenv("USER"),
        "disruption-py",
        ("." if "pytest" in sys.modules else "") + time.strftime("%Y-%m-%d"),
    )
    Path(top).mkdir(parents=True, exist_ok=True)

    # create temporary sub folder
    return mkdtemp(dir=top, prefix=time.strftime("%H.%M.%S-"))


def shot_msg(message: str) -> str:
    """
    Modify a message by prepending a shot format string.

    Parameters
    ----------
    message : str
        The message to modify.

    Returns
    -------
    str
        The modified message to be formatted downstream.
    """
    return "#{shot} | " + message


def shot_msg_patch(mylogger: Type[logger], shot: int):
    """
    Patch a logger by prepending the shot number to its message.

    Parameters
    ----------
    mylogger:
        The logger to modify.
    shot: int
        The shot id to prepend.
    """
    return mylogger.patch(
        lambda r: r.update(message=shot_msg(r["message"]).format(shot=shot))
    )


def get_elapsed_time(elapsed: float) -> str:
    """
    Convert elapsed seconds into human-readable format.

    Parameters
    ----------
    elapsed : float
        Elapsed number of seconds.

    Returns
    -------
    str
        Human-readable formatted message.
    """

    out = []
    d = elapsed // (24 * 3600)
    elapsed -= d * (24 * 3600)
    h = elapsed // 3600
    elapsed -= h * 3600
    m = elapsed // 60
    elapsed -= m * 60
    s = int(elapsed)
    ms = (elapsed - s) * 1000
    if d > 0:
        out += [f"{d:.0f}d"]
    if h > 0:
        out += [f"{h:.0f}h"]
    if m > 0:
        out += [f"{m:.0f}m"]
    if s > 0:
        out += [f"{s:.0f}s"]
    if not out:
        out += [f"{ms:.0f}ms"]
    return " ".join(out)


def to_tuple(
    data: Dict[str, np.ndarray], dim: str
) -> Dict[str, Tuple[str, np.ndarray]]:
    """
    Recreate a dictionary by making all values a 2-tuple with a given string.

    Parameters
    ----------
    data : Dict[str, np.ndarray]
        Dictionary of array data.
    dim : str
        String to be added as first element of the tuple.

    Returns
    -------
    Dict[str, Tuple[str, np.ndarray]]
    """
    return {k: (dim, v) for k, v in data.items()}
