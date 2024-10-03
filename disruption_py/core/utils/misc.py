#!/usr/bin/env python3

"""
Module for utility functions related to class instantiation, data manipulation, and version control.
"""

import subprocess
import warnings
from typing import List

import numpy as np
import pandas as pd


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
