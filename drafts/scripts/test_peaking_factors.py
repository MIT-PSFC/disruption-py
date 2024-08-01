import logging
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data
import matplotlib.pyplot as plt
import MDSplus as mds

shotno = 1150928025
#shotno = 1120621021
signals = ['ip', 'kappa', 'te_peaking', 'te_peaking_ece', 'prad_peaking']

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",	# use the set efit's timebase
    efit_nickname_setting="efit18",	# set the efit
    run_tags=[],
    run_methods=["_get_peaking_factors_ece"],
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

fig, axs = plt.subplots(5, 1, sharex=True)
axs[0].plot(data['time'], data['ip'], label='ip')
axs[1].plot(data['time'], data['kappa'], label='$\kappa$')
axs[2].scatter(te0_time, te0, marker='.', s=1, label='Te0 ECE')
axs[3].scatter(data['time'], data['prad_peaking'], marker='.', label='Prad PF')
axs[4].scatter(data['time'], data['te_peaking_ece'], c='b', marker='o', label='Te PF ECE')
axs[4].scatter(data['time'], data['te_peaking'], c='r', marker='.', label='Te PF TS')
for ax in axs:
    ax.legend()
plt.show()
