#!/usr/bin/env python3

"""
Program to read cmod data for the injection shotlist.
Last Major Update: Henry Wietfeldt (9/09/24)
"""

import logging

import numpy as np
import pandas as pd

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data

FN_IN = 'drafts/ufo_data_in/flattop_disruption_shotlist_2012_to_2016.txt'
FN_OUT = 'drafts/ufo_data_out/cmod_data_2012_to_2016.h5'
DATA_KEY = "sql_data"

shot_list = np.genfromtxt(fname=FN_IN, dtype=int)
print(shot_list.shape)

# default method for pulling disruption-py data
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",	# use the set efit's timebase
    efit_nickname_setting="efit18",	# set the efit
    run_tags=["all"],
    run_columns=[],
    run_methods=[],
    only_requested_columns=False
)

data = get_shots_data(
    tokamak="cmod",
    shotlist_setting=shot_list,
    retrieval_settings=retrieval_settings,
    log_settings=LogSettings(console_log_level=logging.DEBUG),
    output_setting="dataframe",
    num_processes = 20,
)
#data.to_hdf(FN_OUT, key=DATA_KEY, mode='w')