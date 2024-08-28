import logging
import numpy as np
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data
import matplotlib.pyplot as plt
import MDSplus as mds


PLOT = True
#shotno = 1140605022
shotno = 1120828014
#shotno = 1150928025
#shotno = 1140702004 # 1.22 and 1.24 show profile narrowing, 1.2696 shows profile broadening as core cools
#shotno = 1120830026 # 0.52 and 0.57
#shotno = 1150605026 # example of GPC off
#shotno = 1140226024 # example of cutoff
# shotno = 1150605023
#shotno = 1140717002 # Example with LH heating
signals = ['ip', 'kappa', 'p_rad', 'te_peaking', 'te_core_vs_avg_ece', 'te_edge_vs_avg_ece', 'prad_peaking', 'te_width_ece']

# shot_list = np.loadtxt("cmod_shots_2012_to_2016.txt", dtype=int)
# print(shot_list[:5])
# print(shot_list.shape)

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",	# use the set efit's timebase
    efit_nickname_setting="efit18",	# set the efit
    run_tags=[],
    run_methods=["_get_te_profile_params_ece"],
    run_columns=signals,
    only_requested_columns=True
)

data = get_shots_data(
    shotlist_setting=[shotno], # Example shot. Fails for the ~400 shots in my full shotlist
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes = 20,
)

data.to_csv('peaking_factors.csv')
print(data.columns)

if (PLOT):
    # Get Te0 to compare sawteeth to Te peaking
    te0_data = mds.Tree('cmod', shotno).getNode('\\gpc2_te0').getData()
    te0 = te0_data.data()
    te0_time = te0_data.dim_of().data()

    fig, axs = plt.subplots(7, 1, sharex=True)
    axs[0].plot(data['time'], data['ip']/1e6, label='ip', c='k')
    axs[1].plot(te0_time, te0, marker='.', ms=1, label='Te0 ECE', c='k')
    axs[2].plot(data['time'], data['p_rad']/1e6, marker='.', label='Prad [MW]', c='k')
    axs[3].plot(data['time'], data['prad_peaking'], marker='.', label='Prad PF', c='k')
    axs[4].plot(data['time'], data['te_core_vs_avg_ece'], c='b', marker='o', label='ECE')
    axs[4].plot(data['time'], data['te_peaking'], c='r', marker='.', label='TS')
    axs[5].plot(data['time'], data['te_edge_vs_avg_ece'], c='b', marker='o', label='Te Edge')
    axs[6].plot(data['time'], data['te_width_ece'], c='b', marker='o')
    axs[6].set_xlabel('Time [s]', fontsize=14)
    axs[6].set_xlim(0, 1.5)
    axs[6].set_ylim(0, 0.15)
    axs[0].set_ylabel("Ip [MA]")
    axs[1].set_ylabel("Te0 [keV]")
    axs[2].set_ylabel("Prad [MW]")
    axs[3].set_ylabel("Prad PF")
    axs[4].set_ylabel("Te PF")
    axs[5].set_ylabel("Te edge/avg")
    axs[6].set_ylabel("Te width [m]")
    axs[4].legend()
    axs[0].set_title('Cmod Shot ' + str(shotno))
    plt.show()
