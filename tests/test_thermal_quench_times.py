"""
Test automated TQ labels compared to manual labels
Author: Henry Wietfeldt
Source of Manual labels: Henry Wietfeldt
"""

import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

db_man = pd.read_csv('manual_tq_labels.csv')
db_man['shot'] = db_man['shot'].astype(int)
shotlist = db_man['shot'].to_list()
print(f"Number of shots: len(shotlist)")

signals = [
    "ip",
    "thermal_quench_time",
    "t_disrupt"
]

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",  # use the set efit's timebase
    efit_nickname_setting="analysis",  # set the efit
    run_methods=[],
    run_columns=signals,
    only_requested_columns=True,
)

db_auto = get_shots_data(
    shotlist_setting=shotlist,
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_level=logging.WARNING),
    output_setting="dataframe",
    num_processes=os.cpu_count(),
)

# Using pandas because this testing script has not been updated to xarray

# Output test results
db_auto = db_auto.drop_duplicates(subset='shot').drop(columns='time')
db_auto.to_csv('test_thermal_quench_auto_labels.csv')
db_test = db_man.copy()
db_test = db_test.sort_values(by='shot')
db_auto = db_auto.sort_values(by='shot')
db_test = pd.merge(db_test, db_auto, how='outer', on='shot')
db_test = db_test.rename(columns={'thermal_quench_time': 'tq_auto'})
print(db_test)

db_test['within_tq_range'] = (db_test['tq_auto'] > db_test['tq_onset_manual']-0.001) & (db_test['tq_auto'] < db_test['tq_end_manual']+0.001)
db_test['onset_error_s'] = db_test['tq_auto'] - db_test['tq_onset_manual']
db_test = db_test.drop(columns=['notes'])
db_test.to_csv('test_thermal_quench_results.csv')

# Print summary statistics
error = db_test['onset_error_s'].to_numpy()
print(f"Mean |Error| = {1e3*np.mean(np.abs(error)):.3f} ms")
print(f"Median |Error| = {1e3*np.median(np.abs(error)):.3f} ms")
print(f"Std Dev |Error| = {1e3*np.std(np.abs(error)):.3f} ms")
print(f"Min Error = {1e3*np.min(error):.3f} ms")
print(f"Max Error = {1e3*np.max(error):.3f} ms")
print(f"Num Outliers (outisde TQ [start, end] by >1 ms) = {np.sum(~db_test['within_tq_range'])} out of {len(shotlist)} shots")

# Plot onset errors
plt.hist(db_test['onset_error_s'], bins=50)
plt.xlabel('Error in TQ Onset Time (Auto - Manual) [s]', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.show()