#!/usr/bin/env python3

"""This module provides utility functions for pytest like extracting parameters
from pytest command-line arguments and saving data to temporary CSV files."""

import re

import pandas as pd


def extract_param(config):
    """Extract the data column from the pytest command.

    E.g. will return ip given
    `pytest -s tests/test_against_sql.py -k test_data_columns[ip]`

    Params:
        config: pytestconfig fixture

    Returns:
        List[str], the data column if it exists, otherwise None.
    """
    args = config.invocation_params.args
    if len(args) == 0:
        return None
    m = re.search(r"\[(.+)\]$", args[-1])
    param = [m.group(1)] if m is not None else None
    return param


def save_to_csv(data, test_file_path_f, data_source_name):
    """Save a dataframe of fresh or cached data to the tmp testing directory"""
    for shot_id in data:
        pd.DataFrame(data[shot_id]).to_csv(
            test_file_path_f(f"-{data_source_name}-{shot_id}.csv")
        )
