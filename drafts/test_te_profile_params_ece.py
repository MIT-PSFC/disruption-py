"""
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from disruption_py.settings import LogSettings, RetrievalSettings
from disruption_py.workflow import get_shots_data


def read_dpy(shot_list, local_data_dir, efit_name, bypass=False, save_data=True, only_columns=None):
	'''
	saves requested data to dpy_{efit_name}.h5 in the directory specified by local_data_dir
	  - if data file already exists locally, just pulls data from that local file instead
      - can bypass reading existing local data if bypass=True
	  - can skip saving data if save_data = False
	By default grabs all columns ouptut by disruption_py
	  - can grab only certain columns by setting 'only_columns' to a list of requested columns
	returns data as a pandas DataFrame
	
	Adapted from A.R. Saperstein
	'''
	# need to make sure list is of type 'int' and not 'numpy.int64'
	shot_list = np.array(shot_list).tolist()
	
	data_path = f"{local_data_dir}dpy_{efit_name}.h5"	# location of saved data
	data_key = "sql_data"

	# check if requested data already stored locally
	if os.path.exists(data_path) and not bypass:
		data = pd.read_hdf(data_path, key=data_key)
		return data

	# default method for pulling disruption-py data
	retrieval_settings = RetrievalSettings(
		time_setting="disruption_warning",	# use the set efit's timebase
		efit_nickname_setting=efit_name,	# set the efit
		run_methods=None,
		run_columns=only_columns,
		only_requested_columns= only_columns is not None
	)
	data = get_shots_data(
	    shotlist_setting=shot_list,
	    retrieval_settings=retrieval_settings,
	    log_settings=LogSettings(console_level="DEBUG"),
	    output_setting="dataframe",
	 	num_processes = 20,
	)
	# Reset the loguru settings
	if save_data:
		if not os.path.exists(local_data_dir):
			os.mkdir(local_data_dir)
		data.to_hdf(data_path, key=data_key, mode='w')
	return data

if __name__=='__main__':
	
	shot_list = [1160826029]
	time_slice = 1.44
	local_data_dir = '/home/henrycw/projects/disruption-py/drafts/local_data/'
	efit_name='analysis'

	# Pass time slice to be read by physics method
	with open(f'{local_data_dir}time_slice.txt', mode='w') as file:
		file.write(f'{time_slice}')

	df = read_dpy(shot_list=shot_list, local_data_dir=local_data_dir, efit_name=efit_name,
					bypass=True)
	print(df[['time', 'te_width_ece']])
	
	# te_prof = pd.read_csv(f'{local_data_dir}te_prof.csv')

	# fig, axs = plt.subplots(3, 2, sharex=False)
	# axs[0][0].plot(df['time'], df['te_width_ece'], marker='o')
	# axs[0][1].plot(df['time'], df['te_core_vs_avg_ece'], marker='o')
	# axs[0][2].plot(df['time'], df['te_edge_vs_avg_ece'], marker='o')
	# axs[1][0].plot(te_prof['r'], te_prof['te'], marker='o', c='k', linestyle='none')
	# axs[1][0].plot(te_prof['r'], te_prof['te_fit'], c='r')
	# plt.show()

	
	