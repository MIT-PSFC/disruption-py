#!/usr/bin/env python3

"""
Program to plot various quantities used in calculating the thermal quench time 
for a particular shot to compare various methods.
Last Major Update: Henry Wietfeldt (10/04/24)
"""

import logging

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDSplus as mds

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

SHOT_ID = 1120105021
signals = [
    "ip",
    "zcur",
    "wmhd",
    "time_until_disrupt",
    "thermal_quench_time"
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
    shotlist_setting=[SHOT_ID],
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes=20,
)
data.to_csv("thermal_quench_time.csv")
# print(data.columns)

# sxr_data = pd.read_csv('sxr_out.csv')
# fig, axs = plt.subplots(4, 1, sharex=True)
# axs[0].plot(data['time'], np.abs(data['ip']), marker='o')
# axs[1].plot(data['time'], data['zcur'], marker='o')
# axs[2].plot(data['time'], data['wmhd'], marker='o')
# axs[3].scatter(sxr_data['time'], sxr_data['sxr'], marker='.', s=5)
# for ax in axs:
#     ax.axvline(data['thermal_quench_time'][0], linestyle='--', c='r')
#     ax.axvline(data['time'].values[-1], linestyle='--', c='k')
# print("PLOTTING")
# plt.show()

