""" Used in the documentation for the shot data request. """
from typing import List
from disruption_py.settings.shot_data_request import ShotDataRequestParams


def cached_method_params_function(parent_object, shot_data_request_params : ShotDataRequestParams, **kwargs) -> List[str]:
	"""Functions may be passed as some of the parameters of the `parameter_cached_method` and `cached_method` decorators. 
    
    These functions are called at runtime before the decorated method is run. They allow for the dynamic determination of
    parameters such as `used_trees` that may not be defined before runtime.

	Parameters
	----------
	parent_object : ShotDataRequest
		The object that contains the decorated method.
	shot_data_request_params : ShotDataRequestParams
		The ShotDataRequestParams used when calling the decorated method.

	Returns
	-------
	List[str]
		Value of the expected type for the parameter (cuurently all parameters that allow functions are lists of strings)
	"""
	pass