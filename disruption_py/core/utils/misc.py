#!/usr/bin/env python3

"""
Module for utility functions related to class instantiation, data manipulation, and version control.
"""

import os
import subprocess
import time
import warnings
from functools import lru_cache
from pathlib import Path
from tempfile import mkdtemp
from typing import List, Type

import numpy as np
import pandas as pd
from loguru import logger


def instantiate_classes(lst: List):
    """
    Instantiate all classes in a list of classes and objects.

    Parameters
    ----------
    lst : List
        List to instantiate classes from.

    Returns
    -------
    List
        The list with all classes instantiated.
    """
    return [x() for x in lst if isinstance(x, type)]


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


def safe_df_concat(base_df: pd.DataFrame, new_dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Safely concatenate a base DataFrame with a list of new DataFrames.

    Parameters
    ----------
    base_df : pd.DataFrame
        The base DataFrame to concatenate with.
    new_dfs : List[pd.DataFrame]
        A list of new DataFrames to concatenate.

    Returns
    -------
    pd.DataFrame
        The concatenated DataFrame.
    """
    if isinstance(new_dfs, pd.DataFrame):
        new_dfs = [new_dfs]

    all_cols = set(base_df.columns).union(*[set(new_df.columns) for new_df in new_dfs])

    new_dfs = [new_df.dropna(axis=1, how="all") for new_df in new_dfs]
    new_dfs = [
        new_df
        for new_df in new_dfs
        if not new_df.empty and not new_df.isna().all().all()
    ]

    if len(new_dfs) == 0:
        return base_df

    if base_df.empty:
        concat_df = pd.concat(new_dfs, axis=0, ignore_index=True, sort=False)
    else:
        concat_df = pd.concat(
            [base_df] + new_dfs, axis=0, ignore_index=True, sort=False
        )

    missing_cols = all_cols - set(concat_df.columns)
    for col in missing_cols:
        concat_df[col] = np.nan
    return concat_df


@lru_cache
def get_commit_hash() -> str:
    """
    Retrieve the current Git commit hash.

    Returns
    -------
    str
        The short commit hash if available, otherwise None.
    """
    try:
        commit_hash = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .decode("ascii")
            .strip()
        )
    except subprocess.CalledProcessError:
        commit_hash = None
    return commit_hash


@lru_cache
def get_temporary_folder() -> Path:
    """
    Create and return a temporary folder.
    The result is cached to return the same path for different invocations.

    Returns
    -------
    Path
        Resulting temporary folder.
    """

    # create temporary top folder
    top = os.path.join("/tmp", os.getenv("USER"), "disruption-py")
    Path(top).mkdir(parents=True, exist_ok=True)

    # create temporary sub folder
    return mkdtemp(prefix=time.strftime("%Y%m%d-%H%M%S-"), dir=top)


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
