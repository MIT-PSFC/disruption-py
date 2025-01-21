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
import matplotlib.pyplot as plt
import pandas as pd

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

input_fn = 'drafts/scripts/tq_man_labeled_dataset.csv'
output_fn = 'drafts/scripts/test_thermal_quench_onset_output.csv'
physics_method_script = 'disruption_py/machine/cmod/physics.py'
time_above_threshold_scan = [0.001]
normalized_threshold_scan = [0.7]

def modify_script(target_script, time_above_threshold, normalized_threshold):
    with open(target_script, "r") as file:
        lines = file.readlines()
    # Update the parameters in the target script
    for i, line in enumerate(lines):
        stripped_line = line.lstrip()
        if stripped_line.startswith("time_above_threshold ="):
            indent_length = len(line) - len(stripped_line)
            indent = line[:indent_length]
            lines[i] = f"{indent}time_above_threshold = {time_above_threshold}\n"
        elif stripped_line.startswith("normalized_threshold ="):
            indent_length = len(line) - len(stripped_line)
            indent = line[:indent_length]
            lines[i] = f"{indent}normalized_threshold = {normalized_threshold}\n"
    # Write the modified lines back to the target script
    with open(target_script, "w") as file:
        file.writelines(lines)

manual_db = pd.read_csv(input_fn)
np.random.seed(42) # For reproducible output
training_set = np.random.choice(manual_db['shot'].to_numpy(), size=70, replace=False)
output_db = manual_db[manual_db['shot'].isin(training_set)]

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


for tat in time_above_threshold_scan:
    for nt in normalized_threshold_scan:
        modify_script(target_script=physics_method_script, 
                      time_above_threshold=tat,
                      normalized_threshold=nt
        )
        data = get_shots_data(
            shotlist_setting=training_set,
            retrieval_settings=retrieval_settings,
            log_settings=LogSettings(console_log_level=logging.ERROR),
            output_setting="dataframe",
            num_processes=20,
        )
        data = data.sort_values(by=['shot', 'time'])
        data = data.drop_duplicates(subset='shot', keep='first').drop(columns='time').rename(columns={'thermal_quench_time_onset':'tq_onset_auto'})
        data['time_above_threshold'] = tat
        data['normalized_threshold'] = nt
        output_db = output_db.merge(data, how='left', on='shot')
        output_db['tq_onset_auto'] = output_db['tq_onset_auto'].round(5)

output_db.to_csv(output_fn, index=False, columns=['shot','tq_onset_manual','tq_end_manual','tq_onset_auto', 'time_above_threshold',
                                                   'normalized_threshold', 'commit_hash', 'notes'])
