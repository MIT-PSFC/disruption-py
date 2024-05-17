#!/usr/bin/env python3
import logging

import pandas as pd
import disruption_py
from disruption_py.database import create_d3d_handler
from disruption_py.shots import D3DShot

"""
This script demonstrates how to use disruption_py to generate a csv shot dataset
with chosen columns.
"""
FEATURE_COLUMNS = ["time", "ip"]  # Shot columns we want


def generate_full_dataset():
    # Create database handler for grabbing shots from SQL database
    handler = create_d3d_handler()
    # Get all shots from database (you can pass a list of shot_ids to get a subset)
    shots = handler.get_shots()
    # Combine into one dataframe
    df = pd.concat([shot.data[FEATURE_COLUMNS] for shot in shots])
    # Save to csv
    df.to_csv("d3d_shot_data.csv")


def generate_subset_dataset(shot_ids):
    """
    Same as generate_full_dataset but for a subset of shots
    """
    handler = create_d3d_handler()
    shots = handler.get_shots(shot_ids)
    df = pd.concat([shot.data[FEATURE_COLUMNS] for shot in shots])
    df.to_csv("d3d_shot_data.csv")


if __name__ == "__main__":
    logger = logging.getLogger("disruption_py")
    # Output to terminal
    # ch = logging.StreamHandler()

    # Output to file:
    ch = logging.FileHandler("d3d_example.log")

    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    # generate_subset_dataset(['175552','175553'])
    shot_ids = ["191914", "191786"]
    shots = []
    for shot_id in shot_ids:
        shots.append(D3DShot(shot_id, "efit01"))
    df = pd.concat([shot.data[FEATURE_COLUMNS] for shot in shots])
    df.to_csv("d3d_shot_data_local.csv")
