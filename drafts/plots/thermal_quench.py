#!/usr/bin/env python3

"""
Program to plot various quantities used in calculating the thermal quench time 
for a particular shot to compare various methods.
Last Major Update: Henry Wietfeldt (10/04/24)
"""

import logging
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')  
import matplotlib.pyplot as plt
import pandas as pd
import pickle

import disruption_py 
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

#SHOT_ID = 1140515015
#SHOT_ID = 1140827029
SHOT_ID = 1120717002
signals = [
    "ip",
    "zcur",
    "thermal_quench_time_onset"
]

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",  # use the set efit's timebase
    efit_nickname_setting="efit18",  # set the efit
    run_methods=[],
    run_columns=signals,
    only_requested_columns=True,
)

data = get_shots_data(
    shotlist_setting=[SHOT_ID],
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes=1,
)
print(data)

with open('sxr.pkl', 'rb') as f:
    df = pickle.load(f)
df['ip'] = np.abs(df['ip']/1e6)

plt.rcParams['font.size'] = 16
fig, axs = plt.subplots(5, 1, sharex=True, figsize=(10,8))
axs[0].set_xlim(0.2, 0.7)
axs[0].scatter(df['magtime'], df['ip'], marker='.', s=10, c='k')
axs[1].plot(df['magtime'], (df['ip_growth_rate']), marker='.', c='k')
axs[2].scatter(df['t_sxr'], df['core_sxr'], marker='.', s=5, c='k')
axs[3].scatter(df['t_sxr'], df['dcore_sxr_dt'], marker='.', s=5, c='k')
axs[4].scatter(df['efit_time'], df['z0'], marker='o', s=10, c='k')
for ax in axs:
    ax.axvline(df['t_disrupt'], linestyle='--', c='b', label='CQ')
    ax.axvline(df['t_start'], linestyle='--', c='k', label='tstart')
axs[0].set_title('C-Mod Shot: ' + str(SHOT_ID))
axs[0].set_ylabel('Ip [MA]')
axs[1].set_ylabel('dIp/dt [MA/s]')
axs[1].set_ylim(-100,100)
axs[2].set_ylabel('max(SXR) [$W/m^2$]')
axs[3].set_ylabel("dcore_sxr_dt")
axs[3].set_ylim(-8e3, 2e3)
axs[4].set_ylabel('Z0 [m]')
axs[4].set_xlabel("Time [s]")
axs[2].legend(fontsize=12)

plt.show()