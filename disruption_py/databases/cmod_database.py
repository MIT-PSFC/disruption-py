import os
import logging
from disruption_py.databases.database import ShotDatabase
import disruption_py.data
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

CMOD_PROTECTED_COLUMNS = ['dbkey', 'shot', 'time', 'time_until_disrupt', 'ip_error', 'dip_dt', 'beta_p', 'beta_n', 'li', 'n_equal_1_normalized', 'z_error', 'v_z', 'z_times_v_z', 'kappa', 'pressure_peaking', 'H98', 'q0', 'qstar', 'q95', 'dn_dt', 'p_rad_slow', 'p_oh_slow', 'p_icrf', 'p_lh', 'radiated_fraction', 'power_supply_railed', 'v_loop_efit', 'lower_gap', 'upper_gap', 'dbetap_dt', 'dli_dt',
                          'ip', 'zcur', 'n_e', 'dipprog_dt', 'v_loop', 'p_rad', 'p_oh', 'ssep', 'dWmhd_dt', 'dprad_dt', 'Te_width', 'Greenwald_fraction', 'intentional_disruption', 'Te_width_ECE', 'Wmhd', 'n_over_ncrit', 'n_equal_1_mode', 'Mirnov', 'Mirnov_norm_btor', 'Mirnov_norm_bpol', 'Te_peaking', 'ne_peaking', 'Te_peaking_ECE', 'SXR_peaking', 'kappa_area', 'I_efc', 'SXR', 'H_alpha', 'Prad_peaking_CVA', 'commit_hash']

class CModDatabase(ShotDatabase):
	logger = logging.getLogger('disruption_py')

	def __init__(self, driver, driver_file, host, db_name, user, passwd, **kwargs):
		super().__init__(driver, driver_file, host, db_name, user, passwd, protected_columns=CMOD_PROTECTED_COLUMNS, **kwargs)
  
	def default(**kwargs):
		USER = os.getenv('USER')
		# TODO: Catch error if file not found and output helpful error message
		with open(f"/home/{USER}/logbook.sybase_login", "r") as profile:
			content = profile.read().splitlines()[1:]
			db_server = content[0]
			db_name = content[1]
			db_username = content[2]
			#assert db_username == USER, f"db_username:{db_username};user:{USER}"
			db_password = content[3]
		with importlib_resources.path(disruption_py.data, "sqljdbc4.jar") as p:
			db_driver_path = str(p)  # Absolute path to jar file
		CModDatabase.logger.debug(db_driver_path)
		CModDatabase.logger.debug(os.getcwd())
		return CModDatabase(
      		driver="{ODBC Driver 18 for SQL Server}", 
        	driver_file=db_driver_path, 
         	host=db_server,
        	db_name=db_name,
        	user=db_username, 
        	passwd=db_password,
        	**kwargs
        )