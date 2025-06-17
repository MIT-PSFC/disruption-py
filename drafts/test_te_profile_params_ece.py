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
	    log_settings=LogSettings(console_level="WARNING"),
	    output_setting="dataframe",
	 	num_processes = 20,
	)
	# Reset the loguru settings
	if save_data:
		if not os.path.exists(local_data_dir):
			os.mkdir(local_data_dir)
		data.to_hdf(data_path, key=data_key, mode='w')
	return data

def manually_check_profiles():
	
	shot_list = [1100910023]
	time_slices = [0.55, 0.88, 1.0]
	#time_slices = [1.62, 1.64, 1.66]
	#time_slices = [0.26, 0.28, 0.3]
	local_data_dir = '/home/henrycw/projects/disruption-py/drafts/local_data/'
	efit_name='analysis'

	# Pass time slice to be read by physics method
	np.savetxt(f"{local_data_dir}time_slices.txt", time_slices)

	df = read_dpy(shot_list=shot_list, local_data_dir=local_data_dir, efit_name=efit_name,
					bypass=True, save_data=False)
	te_prof = pd.read_csv(f'{local_data_dir}te_prof.csv')

	fig, axs = plt.subplots(3, 2, sharex='col')
	axs[0][0].plot(df['time'], df['te_width_ece'], marker='o')
	axs[0][0].set_ylabel("Te width [m]")
	axs[1][0].plot(df['time'], df['te_core_vs_avg_ece'], marker='o')
	axs[1][0].set_ylabel("Core PF")
	axs[2][0].plot(df['time'], df['te_edge_vs_avg_ece'], marker='o')
	axs[2][0].set_ylabel("Edge PF")
	axs[2][0].set_xlabel("time [s]")
	for i, t in enumerate(time_slices):
		print(i)
		prof = te_prof[te_prof['time']==t]
		axs[i][1].plot(prof['r'], prof['te'], marker='o', c='k', linestyle='none')
		axs[i][1].plot(prof['r'], prof['te_fit'], c='r', label='Fit')
		#axs[i][1].plot(prof['r'], prof['te_guess'], c='g')
		axs[i][1].axvline(prof['r0'].to_numpy()[0], c='gray', label='R0')
		axs[i][1].axvline(prof['r0'].to_numpy()[0] + prof['a'].to_numpy()[0], c='gray', linestyle='--', label='R0 + a')
		axs[i][1].set_ylabel("Te [keV]")
		axs[i][1].set_title(f"t = {t:.3} s")
		axs[i][1].set_ylim(0, np.max(te_prof['te'].to_numpy())+0.5)

		axs[0][0].plot(df[df['time'].round(4)==t]['time'], df[df['time'].round(4)==t]['te_width_ece'], marker='.', c='r')
	axs[0][1].legend()
	axs[2][1].set_xlabel("R [m]")
	plt.suptitle(f"Te Profile Params from ECE\n C-Mod Shot {shot_list[0]}")
	plt.tight_layout()
	plt.show()

def test_batch_run():
	fn_shot_list = "/home/henrycw/projects/ufo-characterization/shotlists/cmod_shots_2012_to_2016.txt"
	shot_list = np.genfromtxt(fn_shot_list, dtype=int)
	local_data_dir = '/home/henrycw/projects/disruption-py/drafts/local_data/'
	efit_name = 'analysis'
	only_columns = ['te_width_ece', 'te_core_vs_avg_ece', 'te_edge_vs_avg_ece', 'te_width', 'kappa']

	df = read_dpy(shot_list=shot_list, local_data_dir=local_data_dir, efit_name=efit_name,
				bypass=False, only_columns=only_columns)
	te_width_ece = df['te_width_ece'].to_numpy()
	print(f"Fraction of NaN slices: {np.sum(np.isnan(te_width_ece))/len(te_width_ece)}")
	print(np.flip(np.sort(te_width_ece)))
	print(np.nanmax(te_width_ece))
	print(np.nanmin(te_width_ece))
	plt.hist(te_width_ece, bins=50)
	plt.yscale('log')
	plt.ylabel('Count')
	plt.xlabel('Te Width [m]')
	plt.title('Te Width from ECE from Analysis Time Slices\non C-Mod 2012-2016')
	plt.show()
	print(df[df['te_width_ece'] > 0.01].sort_values(by='te_width_ece', ascending=True))


if __name__=="__main__":
	test_batch_run()

	
	