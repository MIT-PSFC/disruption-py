#!/usr/bin/env python3

import warnings
from typing import List

import numpy as np
import pandas as pd


def instantiate_classes(l: List):
    """
    Instantiate all classes in a list of classes and objects.

    Parameters
    ----------
    l : List
        List to instantiate classes from.

    Returns
    -------
    List
        The list with all classes instantiated.
    """
    return [x() for x in l if isinstance(x, type)]


def without_duplicates(l: List):
    """
    Get list without duplicates maintaining order.

    Parameters
    ----------
    l : List
        List to get without duplicates.

    Returns
    -------
    List
        The list l with duplicates removed.
    """
    seen = set()
    return [x for x in l if not (x in seen or seen.add(x))]


def safe_cast(array: np.ndarray, dtype, copy=False):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        return array.astype(dtype, copy=copy)


def safe_df_concat(base_df: pd.DataFrame, new_dfs: List[pd.DataFrame]):
    if isinstance(new_dfs, pd.DataFrame):
        new_dfs = [new_dfs]

    new_dfs = [new_df.dropna(axis=1, how="all") for new_df in new_dfs]
    new_dfs = [
        new_df
        for new_df in new_dfs
        if not new_df.empty and not new_df.isna().all().all()
    ]

    if len(new_dfs) == 0:
        return base_df

    if base_df.empty:
        return pd.concat(new_dfs, axis=1, ignore_index=True, sort=False)
    else:
        return pd.concat([base_df] + new_dfs, axis=1, ignore_index=True, sort=False)
