"""
Program to plot various quantities used in calculating the thermal quench time 
for a particular shot to compare various methods.
Author: Henry Wietfeldt
"""

import logging
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt
import pandas as pd
import pickle

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

MAN_LABEL = False
#SHOT_ID = 1140515015 # Hot VDE
#SHOT_ID = 1140827029
#SHOT_ID = 1120717002
#SHOT_ID = 1051206029
#SHOT_ID = 1160714006
#TODO: Shot 1160714006 having issues (low SXR signal). What do we do about ramp-up?
# TODO: Shot 11405522001 has significant SXR spike when plasma hits wall, after main TQ
# TODO: Search for first time at which dSXR/dt is w/in factor of 2 from max?
SHOT_ID = 1120223007  # Doesn't have current spike, not sure if this is a hot VDE
signals = [
    "ip",
    "zcur",
    "thermal_quench_time"
]

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",  # use the set efit's timebase
    efit_nickname_setting="efit21",  # set the efit
    run_methods=[],
    run_columns=signals,
    only_requested_columns=True,
)

data = get_shots_data(
    shotlist_setting=[SHOT_ID],
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes=1,
)
print(data)
print("Got data")

with open('sxr.pkl', 'rb') as f:
    df = pickle.load(f)
df['ip'] = np.abs(df['ip']/1e6)
print(df['cq_onset_time'])

plt.rcParams['font.size'] = 14
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(14,7))
#axs[0].set_xlim(0.6, 0.64)
axs[0].scatter(df['magtime'], df['ip'], marker='.', s=10, c='k')
axs[1].scatter(df['t_sxr'], df['core_sxr_raw'], marker='.', s=5, c='k')
axs[2].scatter(df['t_sxr'], df['core_sxr'], marker='.', s=5, c='k')
axs[3].scatter(df['t_sxr'], df['core_sxr_growth_rate'], marker='.', s=5, c='k')
# axs[4].scatter(df['efit_time'], df['z0'], marker='o', s=10, c='k')
print("Plotting labeled times")
for ax in axs:
    ax.axvline(df['t_disrupt'], linestyle='-', c='k', label='t_disrupt')
    ax.axvline(df['cq_onset_time'], linestyle='--', c='k', label='CQ Onset')
    #ax.axvline(df['t_start'], linestyle='--', c='k', label='tstart')
    if not MAN_LABEL:
        ax.axvline(df['thermal_quench_time_scalar'], linestyle='-', c='r', label='TQ Onset')
        # ax.axvspan(df['cq_onset_time']-0.005, df['cq_onset_time'], alpha=0.15, color='tab:green', label='TQ Midpoint Search Window')
        ax.axvline(df['t_max_sxr_drop'], linestyle='--', c='g', label='TQ Midpoint')
        # for i, t_tq in enumerate(df['thermal_quench_times']):
        #     if i == 0:
        #         ax.axvline(t_tq, linestyle='-', c='r', label='TQ')
        #     else:
        #         ax.axvline(t_tq, linestyle='-', c='r')
        # for i, t_warn in enumerate(df['thermal_quench_warnings']):
        #     if i == 0:
        #         ax.axvline(t_warn, linestyle='--', c='b', label='TQ warn')
        #     else:
        #         ax.axvline(t_warn, linestyle='--', c='b')
axs[0].set_title('C-Mod Shot: ' + str(SHOT_ID))
axs[0].set_ylabel('Ip [MA]')
axs[1].set_ylabel('SXR raw')
axs[2].set_ylabel('SXR filt')
axs[3].set_ylabel(r"$dSXR/dt$ [Hz]")
# axs[3].set_ylim(-8e3, 2e3)
# axs[4].set_ylabel('Z0 [m]')
# axs[4].set_xlabel("Time [s]")
axs[0].legend(fontsize=12)

# fig, axs = plt.subplots(4, 1, sharex=True, figsize=(14,7))
# axs[0].plot(df['magtime'], df['ip'], marker='.', ms=10, c='k')
# axs[1].plot(df['t_sxr'], df['core_sxr_raw'], marker='.', ms=5, c='k')
# axs[2].plot(data['time'], data['p_rad']/1e6, marker='o', ms=10, c='k')
# axs[3].plot(df['efit_time'], df['z0'], marker='o', ms=10, c='k')

# for ax in axs:
#     ax.axvline(df['t_disrupt'], linestyle='--', c='b', label='t_disrupt (DisruptionPy)')
#     ax.axvline(df['thermal_quench_time_scalar'], linestyle='--', c='r', label='TQ Onset (auto)')

# axs[0].set_title('C-Mod Shot: ' + str(SHOT_ID))
# axs[0].set_ylabel('Ip [MA]')
# axs[1].set_ylabel('SXR raw [a.u.]')
# axs[2].set_ylabel('Prad [MW]')
# axs[3].set_ylabel('Z0 [m]')
# axs[3].set_xlabel("Time [s]")
# # axs[-1].set_xlim(0.68, 0.725)
# axs[0].legend()
# plt.show()

plt.show()