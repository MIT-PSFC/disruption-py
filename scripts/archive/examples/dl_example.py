#!/usr/bin/env python3

import pandas as pd

from disruption_py.io.sql import create_d3d_handler

"""
This script demonstrates how to use disruption_py to generate a csv shot dataset
with chosen columns.
"""
FEATURE_COLUMNS = ["time", "n_e", "ip", "g_f", "z_eff", "bt0"]  # Shot columns we want


def generate_full_dataset():
    # Create database handler for grabbing shots from SQL database
    handler = create_d3d_handler()
    # Get all shots from database (you can pass a list of shot_ids to get a subset)
    shots = handler.get_shots(efit_tree="EFIT01")
    # Combine into one dataframe
    df = pd.concat([shot.data[FEATURE_COLUMNS] for shot in shots])
    # Save to csv
    df.to_csv("dl_example_data.csv")


# Same as generate_full_dataset but for a subset of shots
def generate_subset_dataset(shot_ids):
    handler = create_d3d_handler()
    shots = handler.get_shots(shot_ids, "EFIT01")
    df = pd.concat([shot.data[FEATURE_COLUMNS] for shot in shots])
    df.to_csv("dl_example_data.csv")


if __name__ == "__main__":
    generate_subset_dataset(["175552"])  # 191914
