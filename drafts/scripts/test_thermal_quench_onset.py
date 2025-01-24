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
import MDSplus as mds

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

input_fn = 'drafts/scripts/tq_man_labeled_dataset_large.csv'
output_fn = 'drafts/scripts/test_thermal_quench_onset_large_output.csv'

manual_db = pd.read_csv(input_fn)

#shotlist = manual_db['shot'].to_numpy()
np.random.seed(42) # For reproducible output
training_set = np.random.choice(manual_db['shot'].to_numpy(), size=70, replace=False)
testing_set = np.setdiff1d(manual_db['shot'].to_numpy(), training_set)

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

data = get_shots_data(
    shotlist_setting=training_set,
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes=20,
)

data = data.sort_values(by=['shot', 'time'])
data = data.drop_duplicates(subset='shot', keep='first').drop(columns='time').rename(columns={'thermal_quench_time_onset':'tq_onset_auto'})
output_db = manual_db.merge(data, how='left', on='shot')
output_db['tq_onset_auto'] = output_db['tq_onset_auto'].round(5)
output_db.to_csv(output_fn, index=False, columns=['shot','tq_onset_manual','tq_end_manual','tq_onset_auto', 'commit_hash', 'notes'])
