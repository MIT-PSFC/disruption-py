import logging
from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

shotno = 1160805001
signals = ['ip', 'te_width', 'te_width_ece']

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",	# use the set efit's timebase
    efit_nickname_setting="efit18",	# set the efit
    run_tags=[],
    run_methods=["_get_peaking_factors_ece"],
    run_columns=signals,
    only_requested_columns=False
)

data = get_shots_data(
    shotlist_setting=[shotno], # Example shot. Fails for the ~400 shots in my full shotlist
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes = 20,
)

data.to_csv('peaking_factors.csv')