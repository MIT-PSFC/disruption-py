import os
import logging
from disruption_py.databases import ShotDatabase

D3D_PROTECTED_COLUMNS = []

class D3DDatabase(ShotDatabase):
	logger = logging.getLogger('disruption_py')

	def __init__(self, driver, host, db_name, user, passwd, **kwargs):
		super().__init__(driver, host, db_name, user, passwd, protected_columns=D3D_PROTECTED_COLUMNS, **kwargs)
  
	def default(**kwargs):
		USER = os.getenv('USER')
		# TODO: Catch error if file not found and output helpful error message
		with open(f"/home/{USER}/D3DRDB.sybase_login", "r") as profile:
			content = profile.read().splitlines()[1:]
			db_username = content[0]
			assert db_username == USER, f"db_username:{db_username};user:{USER}"
			db_password = content[1]
		return D3DDatabase(
      		driver="{ODBC Driver 18 for SQL Server}", 
        	driver_file="UNKNOWN", # TODO
         	host= "UNKNOWN" + "database=code_rundb", # TODO
        	db_name= "UNKNOWN", # TODO
        	user=db_username, 
        	passwd=db_password,
        	**kwargs
        )

	def get_efit_tree(self, shot_id):
		with self.conn.cursor() as curs:
			curs.execute(
				f"select tree from plasmas where shot = {shot_id} and runtag = 'DIS' and deleted = 0 order by idx")
			efit_trees = curs.fetchall()
		if len(efit_trees) == 0:
			efit_trees = [('EFIT01',)]
			# with self.conn.cursor() as curs:
				# curs.execute(f"select tree from plasmas where shot = {shot_id} and deleted = 0 order by idx")
				# efit_trees = curs.fetchall()
		efit_tree = efit_trees[-1][0]
		return efit_tree
