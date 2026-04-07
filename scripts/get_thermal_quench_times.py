"""
This script gets the times of the thermal quench onset for a specified shotlist.
Useful for testing the thermal quench labler.
We can delete this when we merge the thermal quench labeler into dev
Author: Henry Wietfeldt
"""

import os
import logging

import numpy as np

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data


# Shotlist of all C-Mod shots 2012-2016
SHOTLIST_FN = '/home/henrycw/projects/ufo-characterization/shotlists/cmod_shots_2012_to_2016.txt'
NUM_SUBSET = 500 # Number of shots to randomly select from Shotlist

shotlist = np.genfromtxt(SHOTLIST_FN, dtype=int)
rng = np.random.default_rng(seed=42)
rng.shuffle(shotlist)
if len(shotlist) <= NUM_SUBSET:
    shots_to_use = shotlist
else:
    shots_to_use = shotlist[:NUM_SUBSET]

# shots_to_use = [1140821020, 1140520016 ,1150710007]

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",  # use the set efit's timebase
    efit_nickname_setting="analysis",  # set the efit
    run_methods=[],
    run_columns=["ip", "thermal_quench_time"],
    only_requested_columns=True,
)

results = get_shots_data(
    shotlist_setting=shots_to_use,
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_level=logging.WARNING),
    num_processes=os.cpu_count(),
)

# Write contents to csv for easy inspection
results.to_dataframe().to_csv('tq_labels.csv')