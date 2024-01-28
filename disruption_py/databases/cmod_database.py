import os
import logging
from disruption_py.databases.database import ShotDatabase
from disruption_py.utils.constants import CMOD_PROTECTED_COLUMNS

class CModDatabase(ShotDatabase):
	logger = logging.getLogger('disruption_py')

	def __init__(self, driver, host, db_name, user, passwd, **kwargs):
		super().__init__(driver, host, db_name, user, passwd, protected_columns=CMOD_PROTECTED_COLUMNS, **kwargs)
  
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
		return CModDatabase(
      		driver="{ODBC Driver 18 for SQL Server}", 
         	host=db_server,
        	db_name=db_name,
        	user=db_username, 
        	passwd=db_password,
        	**kwargs
        )