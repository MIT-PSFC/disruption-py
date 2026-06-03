"""
Plot timetraces of shots showing different disruption times
for testing thermal quench time labeler
"""

import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

USE_PREV_DPY_RESULTS = True
N_SHOTS = 30

df_man = pd.read_csv('manual_tq_labels.csv')
df_man['shot'] = df_man['shot'].astype(int)
shotlist = df_man['shot'].to_list()
print(f"Number of shots: {len(shotlist)}")

signals = [
    "ip",
    "thermal_quench_time",
    "t_disrupt",
    "core_sxr"
]

if USE_PREV_DPY_RESULTS:
    df = pd.read_csv('tq_df.csv')
else:
    # default method for pulling disruption-py data
    retrieval_settings = RetrievalSettings(
        time_setting="disruption_warning",  # use the set efit's timebase
        efit_nickname_setting="analysis",  # set the efit
        run_methods=[],
        run_columns=signals,
        only_requested_columns=True,
    )

    df = get_shots_data(
        shotlist_setting=shotlist,
        retrieval_settings=retrieval_settings,
        log_settings=LogSettings(console_level=logging.WARNING),
        output_setting="dataframe",
        num_processes=os.cpu_count(),
    )
    df.to_csv('tq_df.csv')

df = df.merge(df_man[['shot', 'tq_onset_manual']], on='shot', how='left')

# Select 30 random shots
rng = np.random.default_rng(seed=42)
rng.shuffle(shotlist)
shot_subset = shotlist[:N_SHOTS]

fig, axs = plt.subplots(N_SHOTS, 1, figsize=(5, 14))

for i, s in enumerate(shot_subset):
    df_s = df[df['shot']==s]
    axs[i].plot(df_s['time'], df_s['ip'], c='g')
    axs[i].set_yticks([])
    axs[i].set_xticks([])
    axs[i].set_xlim(0, 2)

plt.show()