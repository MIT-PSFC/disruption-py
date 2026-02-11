#!/usr/bin/env python3

"""
Module for utility functions related to class instantiation, data manipulation, and version control.
"""

import importlib.metadata
import os
import subprocess
import sys
from datetime import datetime
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
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
            )
            .decode("ascii")
            .strip()
        )
    except subprocess.CalledProcessError:
        commit_hash = ""
    return commit_hash


@lru_cache
def get_metadata() -> Dict[str, str]:
    """
    Gather workflow metadata.

    Returns
    -------
    Dict[str, str]
        The workflow metadata.
    """

    package, *_ = __name__.split(".")
    version = importlib.metadata.version(package)
    tag = "v" + ".".join(version.split(".")[:2])
    repo = "https://github.com/MIT-PSFC/disruption-py"
    commit = get_commit_hash()
    if commit:
        source = f"{repo}/tree/{commit}"
    else:
        source = f"{repo}/releases/tag/{tag}"

    metadata = {
        "user": os.getenv("USER"),
        "host": os.uname().nodename,
        "time": datetime.now().isoformat(),
        "package": package,
        "version": version,
        "commit": commit,
        "source": source,
    }
    if not commit:
        metadata.pop("commit")
    return metadata


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
        os.getenv("LOCALSCRATCH", "/tmp"),
        os.getenv("USER"),
        "disruption-py",
        ("." if "pytest" in sys.modules else "") + datetime.now().strftime("%Y-%m-%d"),
    )
    Path(top).mkdir(parents=True, exist_ok=True)

    # create temporary sub folder
    return mkdtemp(dir=top, prefix=datetime.now().strftime("%H.%M.%S-"))


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
