import os
import jaydebeapi
import logging
from disruption_py.utils.math_utils import save_open_plots
from pandas.api.types import is_numeric_dtype
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from disruption_py.shots import D3DShot
import disruption_py.data
from disruption_py.databases import ShotDatabase
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

D3D_PROTECTED_COLUMNS = []

class D3DDatabase(ShotDatabase):
	def __init__(self, driver, driver_file, host, user, passwd, **kwargs):
			super().__init__(driver, driver_file, host +
							"database=d3drdb", user, passwd, protected_columns=D3D_PROTECTED_COLUMNS, **kwargs)
			self.tree_conn = jaydebeapi.connect(self.driver,
												host + "database=code_rundb", 
            									[self.user, self.passwd],
												self.driver_file)
	def default(**kwargs):
		USER = os.getenv('USER')
		# TODO: Catch error if file not found and output helpful error message
		with open(f"/home/{USER}/D3DRDB.sybase_login", "r") as profile:
			content = profile.read().splitlines()
			db_username = content[0]
			assert db_username == USER, f"db_username:{db_username};user:{USER}"
			db_password = content[1]
		with importlib_resources.path(disruption_py.data, "sqljdbc4.jar") as p:
			db_driver_path = str(p)
		return D3DDatabase("com.microsoft.sqlserver.jdbc.SQLServerDriver", db_driver_path, "jdbc:sqlserver://d3drdb.gat.com:8001;", db_username, db_password, **kwargs)
 
	def get_efit_tree(self, shot_id):
		with self.tree_conn.cursor() as curs:
			curs.execute(
				f"select tree from plasmas where shot = {shot_id} and runtag = 'DIS' and deleted = 0 order by idx")
			efit_trees = curs.fetchall()
		if len(efit_trees) == 0:
			efit_trees = [('EFIT01',)]
			# with self.tree_conn.cursor() as curs:
				# curs.execute(f"select tree from plasmas where shot = {shot_id} and deleted = 0 order by idx")
				# efit_trees = curs.fetchall()
		efit_tree = efit_trees[-1][0]
		return efit_tree

	def validate_shot(self, shot_id, visualize_differences=False, output_dir='./'):
		"""
		Compare shot data currently in disruption database to what is calculated by the shot object.

		Parameters
		----------
		shot_id : int or str
			Shot number to validate
		visualize_differences : bool, optional
			Whether to plot the differences between the two dataframes, by default False

		Returns
		-------
		bool
			Whether the shot data is correct according to the disruption database
		"""
		shot_id = int(shot_id)
		try:
			true_shot = self.get_shot(shot_id)
			if true_shot is None:
				logging.debug("Shot not in database")
				return False, "Shot not in database"
		except Exception as e:
			logging.warning(f"Failed to load shot {shot_id} from SQL database")
			logging.debug(e)
			return False, "Failed to load shot"
		local_shot = D3DShot(shot_id, self.get_efit_tree(
			shot_id), disruption_time=self.get_disruption_time(shot_id))
		comparison = pd.DataFrame()
		for col in list(local_shot.data.columns):
			if col in list(true_shot.data.columns):
				if is_numeric_dtype(local_shot.data[col]):
					comparison[col] = np.abs((true_shot.data[col] - local_shot.data[col])/true_shot.data[col])
		validation = False
		if not comparison.empty:
			logging.debug("Shot data is not correct")
			for col in comparison.columns:
				if (comparison[col] > 1e-4).any():
					fig, axs = plt.subplots(2, 1, sharex=True)
					axs[0].plot(
						true_shot.data['time'], true_shot.data[col], label='database')
					axs[0].plot(
						local_shot.data['time'], local_shot.data[col], label='local', linestyle="dotted")
					axs[1].plot(local_shot.data['time'], comparison[col])
					axs[0].set_title(f"Shot {shot_id}:{col}")
					axs[1].set_title("Normalized Error")
					axs[0].legend()
			save_open_plots(output_dir + f"shot_{shot_id}_error_graphs.pdf")
			plt.close('all')
			if visualize_differences:
				plt.show()
		return validation, comparison