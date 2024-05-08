'''
Be sure to run this script using 'disruption-python' instead of 'python' to use
the correct environment.
'''

from disruption_py.handlers import CModHandler
from disruption_py.settings import ShotSettings, LogSettings
import logging
import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd
import os

# load all possile feature-sets and shot-lists
with open("feature_sets.json","r") as f:
    feature_sets = json.load(f)

with open("shot_lists.json","r") as f:
    shot_lists = json.load(f)

shot_list_name = 'random100'
feature_set_name = 'mhd'

# choose feature-set and shot-list
shot_list = np.array(shot_lists[shot_list_name])
features = feature_sets[feature_set_name]

#print(shot_list)

# pull data

if os.path.exists(f"{shot_list_name}-{feature_set_name}.csv"):
    shot_data = pd.read_csv(f"{shot_list_name}-{feature_set_name}.csv")
elif os.path.exists(f"{shot_list_name}-{feature_set_name}.h5"):
    shot_data = pd.read_hdf(f"{shot_list_name}-{feature_set_name}.hdf")
else:
    handler = CModHandler(mds_connection_str='alcdata-archives')   # {'alcdata-archives': MDSplus, 'DoNotConnect': hsds}

    shot_settings = ShotSettings(
        # logging
        log_settings=LogSettings(
            log_file_path=None,
            file_log_level=logging.WARNING,
            log_file_write_mode="w",
            log_to_console=True,
            console_log_level=logging.WARNING,
            use_custom_logging=False,
        ),

        # data settings
        existing_data_request=None,
        efit_tree_name="efit18",
        
        # method selection
        run_methods=[],
        run_tags=[],    
        run_columns=features,   # add list of features here
        only_requested_columns=True,
        shot_data_requests=[],
        
        # timebase settings
        set_times_request = "efit", # use efit timebase
        signal_domain = "full",
        use_existing_data_timebase = False,
        interpolation_method= "linear",
    )

    shot_data = handler.get_shots_data(
        shot_ids_request=shot_list,
        shot_settings=shot_settings,
        output_type_request = "dataframe", # output a list of dataframes
        num_processes=20,
    )

    shot_data.to_csv(f"{shot_list_name}-{feature_set_name}.csv")
    shot_data.to_hdf(f"{shot_list_name}-{feature_set_name}.h5", key='data')

print(shot_data)

# check csv and h5 have the same data
csv_data = pd.read_csv(f"{shot_list_name}-{feature_set_name}.csv")
h5_data = pd.read_hdf(f"{shot_list_name}-{feature_set_name}.h5")

diff = csv_data - h5_data

# check data
for f in features:

    # check if feature was actually pulled
    if f not in shot_data.keys():
        print(f"'{f}' missing from data")
        continue

    # flag features that are just nans
    if len(np.where(~np.isnan(shot_data[f]))[0]) == 0:
        print(f"'{f}' only contrains NaNs!!!")
        continue

    if len(features) < 20:
        fig,ax = plt.subplots()
        bin_num = 100
        bins = 100
        ax.hist(shot_data[f],bins=bins)
        ax.set_yscale('log')
        ax.set_ylabel('counts')
        ax.set_title(f)
    else:
        print(shot_data[f])

print(f"{len(np.unique(shot_data['shot']))} shots in the dataset")

if len(features) < 20:
    plt.show()