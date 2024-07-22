#!/usr/bin/env python3

import pandas as pd


def save_to_csv(data, module_file_path_f, data_source_name):
    """Save a dataframe of MDSPlus or SQL data to the tmp testing directory"""
    for shot_id in data:
        pd.DataFrame(data[shot_id]).to_csv(
            module_file_path_f(f"-{data_source_name}-{shot_id}.csv")
        )
