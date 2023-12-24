# Shot Id Requests
A module for handling shot id requests passed in the [get_shots_data][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] 
function. Shot id requests are used by disruption_py to get the shot ids of shots that should have data retrieved from MDSplus.

This module defines the abstract class `ShotIdRequest` that can have subclasses passed as the shot_id_request
argument to the [get_shots_data][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] function.
It also provides built_in classes and mappings to easily define shot ids for common use cases.

Currently, these are the options that can be passed as the shot_id_request argument to 
[get_shots_data][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data]:

- A subclass of `ShotIdRequest` from the `disruption_py.settings.shot_id_requests` module.
- A single shot id
- A python list of any other shot id request option listed here
- A dictionary key as a string from the built-in mappings to data files in the `_get_shot_id_request_mappings` dictionary: 
	```python
	--8<--
	disruption_py/settings/shot_id_requests.py:get_shot_id_request_dict
	--8<--
	```
- A file path as a string with its suffix mapped to a ShotIdRequest type in the `_file_suffix_to_shot_id_request` dictionary:
	```python
	--8<--
	disruption_py/settings/shot_id_requests.py:file_suffix_to_shot_id_request_dict
	--8<--
	```
- A dictionary mapping tokamak type strings to the desired shot id request option for that tokamak.  E.g. `{'cmod': 'paper'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"

::: disruption_py.settings.shot_id_requests
    handler: python