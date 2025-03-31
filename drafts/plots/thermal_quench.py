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
#SHOT_ID = 1120717002
SHOT_ID = 1140827029
#TODO: Shot 1160714006 having issues (low SXR signal). What do we do about ramp-up?
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

plt.rcParams['font.size'] = 14
fig, axs = plt.subplots(6, 1, sharex=True, figsize=(12,8))
#axs[0].set_xlim(0.6, 0.64)
axs[0].scatter(df['magtime'], df['ip'], marker='.', s=10, c='k')
axs[1].plot(df['magtime'], (df['ip_growth_rate']), marker='.', c='k')
axs[2].scatter(df['t_sxr'], df['core_sxr_raw'], marker='.', s=5, c='k')
axs[3].scatter(df['t_sxr'], df['core_sxr'], marker='.', s=5, c='k')
axs[4].scatter(df['t_sxr'], df['core_sxr_growth_rate'], marker='.', s=5, c='k')
axs[5].scatter(df['efit_time'], df['z0'], marker='o', s=10, c='k')
for ax in axs:
    ax.axvline(df['t_disrupt'], linestyle='--', c='b', label='CQ')
    #ax.axvline(df['t_start'], linestyle='--', c='k', label='tstart')
    for i, t_tq in enumerate(df['thermal_quench_times']):
        if i == 0:
            ax.axvline(t_tq, linestyle='--', c='r', label='TQ')
        else:
            ax.axvline(t_tq, linestyle='--', c='r')
axs[0].set_title('C-Mod Shot: ' + str(SHOT_ID))
axs[0].set_ylabel('Ip [MA]')
axs[1].set_ylabel(r'$\gamma_{ip}$ [Hz]')
axs[1].set_ylim(-300,300)
axs[2].set_ylabel('SXR raw [a.u.]')
axs[2].set_ylabel('SXR [a.u.]')
axs[4].set_ylabel(r"$\gamma_{SXR}$ [Hz]")
axs[4].set_ylim(-8e3, 2e3)
axs[5].set_ylabel('Z0 [m]')
axs[5].set_xlabel("Time [s]")
axs[0].legend(fontsize=12)

plt.show()