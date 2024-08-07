import logging
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data
import matplotlib.pyplot as plt
import MDSplus as mds

#shotno = 1150928025
#shotno = 1120828014
shotno = 1120830026
signals = ['ip', 'kappa', 'p_rad', 'te_peaking', 'te_peaking_ece', 'prad_peaking', 'te_width_ece']

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

# Get Te0 to compare sawteeth to Te peaking
te0_data = mds.Tree('cmod', shotno).getNode('\\gpc2_te0').getData()
te0 = te0_data.data()
te0_time = te0_data.dim_of().data()

fig, axs = plt.subplots(6, 1, sharex=True)
axs[0].plot(data['time'], data['ip']/1e6, label='ip', c='k')
axs[1].plot(te0_time, te0, marker='.', ms=1, label='Te0 ECE', c='k')
axs[2].plot(data['time'], data['p_rad']/1e6, marker='.', label='Prad [MW]', c='k')
axs[3].plot(data['time'], data['prad_peaking'], marker='.', label='Prad PF', c='k')
axs[4].plot(data['time'], data['te_peaking_ece'], c='b', marker='o', label='ECE')
axs[4].plot(data['time'], data['te_peaking'], c='r', marker='.', label='TS')
axs[5].plot(data['time'], data['te_width_ece'], c='b', marker='.')
axs[5].set_xlabel('Time [s]', fontsize=14)
axs[5].set_xlim(0, 1.5)
axs[0].set_ylabel("Ip [MA]")
axs[1].set_ylabel("Te0 [keV]")
axs[2].set_ylabel("Prad [MW]")
axs[3].set_ylabel("Prad PF")
axs[4].set_ylabel("Te PF")
axs[5].set_ylabel("Te width [m]")
axs[4].legend()
fig.suptitle('Cmod Shot ' + str(shotno))
plt.show()
