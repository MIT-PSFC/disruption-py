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
seed=42
training_frac = 0.7

manual_db = pd.read_csv(input_fn)

#shotlist = manual_db['shot'].to_numpy()
manual_db = pd.read_csv(input_fn)
training_set_05 = manual_db[manual_db['shot'] < 1060000000].sample(frac=training_frac, random_state=seed)
training_set_12_to_16 = manual_db[manual_db['shot'] > 1120000000].sample(frac=training_frac, random_state=seed)
training_set = pd.concat([training_set_05, training_set_12_to_16])

# Get testing set
df_helper = manual_db.merge(training_set, how='left', indicator=True)
testing_set = df_helper[df_helper['_merge'] == 'left_only'].drop(columns=['_merge'])

signals = [
    "thermal_quench_time_onset",
    "time_until_disrupt"
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
    shotlist_setting=testing_set['shot'].to_numpy(),
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes=20,
)

data['t_cq'] = data['time'] + data['time_until_disrupt']
data = data.sort_values(by=['shot', 'time'])
data = data.drop_duplicates(subset='shot', keep='first').drop(columns=['time', 'time_until_disrupt']).rename(columns={'thermal_quench_time_onset':'tq_onset_auto'})
print(data)
output_db = manual_db.merge(data, how='outer', on='shot')
output_db['tq_onset_auto'] = output_db['tq_onset_auto'].round(5)
output_db.to_csv(output_fn, index=False, columns=['shot','tq_onset_manual','tq_end_manual','tq_onset_auto', 't_cq', 'commit_hash', 'notes'])
