#!/usr/bin/env python3

"""
Program to test the C-Mod automated off-line thermal quench onset labeling method
(get_thermal_quench_time_onset() physics method) on a manually labeled database.
Database:
Manually labeled using SXR, Te0 from ECE, Wmhd, among other signals by Henry Wietfeldt
Last Major Update: Henry Wietfeldt (01/06/25)
"""

import logging

import numpy as np
import pandas as pd
import yaml
import MDSplus.mdsExceptions as mdse

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

input_fn = 'drafts/scripts/tq_man_labeled_dataset_large.csv'
output_fn = 'drafts/scripts/train_thermal_quench_onset_output.csv'
physics_method_script = 'disruption_py/machine/cmod/physics.py'
min_time_above_threshold_scan = [0.004]
normalized_threshold_scan = [0.5]

def modify_param_file(min_time_above_threshold, normalized_threshold, param_file='tq_params.yaml'):
    """
    Updates the YAML params file with the new threshold parameters.
    """
    # Load current thermal quench params
    with open(param_file, "r") as f:
        tq_params = yaml.safe_load(f)

    # Modify the config
    tq_params["min_time_above_threshold"] = min_time_above_threshold
    tq_params["normalized_threshold"] = normalized_threshold

    # Write it back to disk
    with open(param_file, "w") as f:
        yaml.safe_dump(tq_params, f)

manual_db = pd.read_csv(input_fn)
np.random.seed(42) # For reproducible output
training_set = np.random.choice(manual_db['shot'].to_numpy(), size=70, replace=False)
manual_db = manual_db[manual_db['shot'].isin(training_set)]

signals = [
    "thermal_quench_time_onset"
]

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",  # use the set efit's timebase
    efit_nickname_setting="efit18",  # set the efit
    run_tags=[],
    run_methods=[],
    run_columns=signals,
    only_requested_columns=True,
)

scan_size = len(min_time_above_threshold_scan)*len(normalized_threshold_scan)
trials = pd.MultiIndex.from_product([min_time_above_threshold_scan, normalized_threshold_scan], 
                                    names=['min_time_above_threshold', 'normalized_threshold'])
output_db = pd.DataFrame({'square_loss': np.zeros(scan_size),
                          'square_loss_norm': np.zeros(scan_size),
                          'late_square_loss': np.zeros(scan_size),
                          'outliers': np.zeros(scan_size)
                          }, index=trials)
print(output_db)
for mtat in min_time_above_threshold_scan:
    for nt in normalized_threshold_scan:
        modify_param_file(min_time_above_threshold=mtat, normalized_threshold=nt)
        try:
            data = get_shots_data(
                shotlist_setting=training_set,
                retrieval_settings=retrieval_settings,
                log_settings=LogSettings(console_log_level=logging.WARNING),
                output_setting="dataframe",
                num_processes=20,
            )
        except mdse.MdsException as e:
            print("Error " + str(mtat) + ", " + str(nt))
            print(e)
            continue
        data = data.sort_values(by=['shot', 'time'])
        data = data.drop_duplicates(subset='shot', keep='first').drop(columns='time').rename(columns={'thermal_quench_time_onset':'tq_onset_auto'})
        data['min_time_above_threshold'] = mtat
        data['normalized_threshold'] = nt
        # Types of losses:
        # Square loss, square loss normalized by width, late square loss, outliers (> 1 ms)
        data = manual_db.merge(data, how='left', on='shot')
        data['square_loss'] = (1000*(data['tq_onset_auto'] - data['tq_onset_manual']))**2 # in ms
        data['square_loss_norm'] = data['square_loss'] / (1000*(data['tq_end_manual']-data['tq_onset_manual']))**2
        data['late_square_loss'] = data['square_loss']
        data.loc[data['tq_onset_auto'] < data['tq_onset_manual'], 'late_square_loss'] = 0
        data['outliers'] = 0
        data.loc[np.abs(data['tq_onset_auto'] - data['tq_onset_manual'])/(data['tq_end_manual']-data['tq_onset_manual']) > 2, 'outliers'] = 1
        print(data[['shot', 'tq_onset_manual', 'tq_end_manual', 'tq_onset_auto', 'square_loss', 'square_loss_norm', 'late_square_loss', 'outliers']])
        # Output stats
        output_db.loc[(mtat, nt), 'square_loss'] = data['square_loss'].sum() / len(training_set)
        output_db.loc[(mtat, nt), 'square_loss_norm'] = data['square_loss_norm'].sum() / len(training_set)
        output_db.loc[(mtat, nt), 'late_square_loss'] = data['late_square_loss'].sum() / len(training_set)
        output_db.loc[(mtat, nt), 'outliers'] = data['outliers'].sum() / len(training_set)
        output_db.to_csv(output_fn, index=True)
        print(output_db)
