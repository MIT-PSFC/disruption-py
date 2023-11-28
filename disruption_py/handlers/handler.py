from enum import Enum
from disruption_py.utils.mappings.tokemak_helpers import get_shot_class_for_shot_id, get_database_for_shot_id

class Handler:
	
	@staticmethod
	def get_sql_data_for_shot(shot_id, **kwargs):
		shot_database_class = get_database_for_shot_id(shot_id)
		shot_database = shot_database_class.default(**kwargs)
		return shot_database.get_shot_data(shot_ids=[shot_id])

	@staticmethod
	def get_mdsplus_data_for_shot(shot_id, **kwargs):
		shot_class = get_shot_class_for_shot_id(shot_id)
		shot = shot_class(shot_id=shot_id, **kwargs)
		return shot.data