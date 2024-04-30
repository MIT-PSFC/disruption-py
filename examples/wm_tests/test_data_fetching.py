'''
Be sure to run this script using 'disruption-python' instead of 'python' to use
the correct environment.
'''

from disruption_py.databases.cmod_database import CModDatabase
from disruption_py.handlers import CModHandler
from disruption_py.settings import ShotSettings, LogSettings
import logging
import numpy as np
import time
import matplotlib.pyplot as plt
import pickle

# load all possile feature-sets and shot-lists
# - try xarray instead
with open('shot_lists.pkl', 'rb') as f:
	shot_lists = pickle.load(f)
with open('feature_sets.pkl', 'rb') as f:
	feature_sets = pickle.load(f)

# choose feature-set and shot-list
shot_list = np.array(shot_lists['random100'])
features = feature_sets['mhd']

# pull data

start_time = time.time()

handler = CModHandler()

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
    efit_tree_name="analysis",
    
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

print()
print(shot_data.keys())

end_time = time.time()

print(f"Elapsed time: {end_time-start_time} seconds")

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