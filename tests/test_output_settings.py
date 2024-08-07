#!/usr/bin/env python3

import os

import pandas as pd

from disruption_py.workflow import get_shots_data


def test_output_exists(shotlist, test_file_path_f):
    """Test creation of all output formats."""
    results = get_shots_data(
        shotlist_setting=shotlist,
        output_setting=[
            "list",
            "dataframe",
            "dict",
            test_file_path_f(".csv"),
            test_file_path_f(".hdf5"),
        ],
        num_processes=2,
    )

    list_output, df_output, dict_output, csv_processed, hdf_processed = results
    assert isinstance(list_output, list)
    assert isinstance(df_output, pd.DataFrame)
    assert isinstance(dict_output, dict)
    assert csv_processed == hdf_processed == len(list_output) == len(shotlist)
    assert os.path.exists(test_file_path_f(".csv"))
    assert os.path.exists(test_file_path_f(".hdf5"))
