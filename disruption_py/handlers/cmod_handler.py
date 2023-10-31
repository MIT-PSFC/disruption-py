
from disruption_py.handlers.handler import Handler
from disruption_py.databases import CModDatabase
from disruption_py.shots import CModShot

class CModHandler:
	
	def __init__(self, database : CModDatabase = None, **kwargs):
		if database is None:
			self.database : CModDatabase = CModDatabase.default()
		else:
			self.database : CModDatabase = database
  
	def get_shot(self, shot_id, **kwargs):
		"""
		Get shot data from CMOD.
		"""
		shot_data = self.database.get_shot_data(shot_ids=[shot_id])
		disruption_time=self.database.get_disruption_time(shot_id)
		return CModShot(shot_id=shot_id, existing_data=shot_data, disruption_time=disruption_time, **kwargs)