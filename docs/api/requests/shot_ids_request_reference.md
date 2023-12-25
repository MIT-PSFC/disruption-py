A module for handling shot ids requests passed in the [get_shots_data][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] 
function. Shot ids requests are used by disruption_py to get the shot ids of shots that should have data retrieved from MDSplus.

This module defines the abstract class `ShotIdsRequest` that can have subclasses passed as the shot_ids_request
argument to the [get_shots_data][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] function.
It also provides built_in classes and mappings to easily define shot ids for common use cases.

Currently, these are the options that can be passed as the shot_ids_request argument to 
[get_shots_data][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data]:

- A subclass of `ShotIdsRequest` from the `disruption_py.settings.shot_ids_requests` module.
- A single shot id
- A python list of any other shot id request option listed here
- A dictionary key as a string from the built-in mappings to data files in the `_get_shot_ids_request_mappings` dictionary: 
	```python
	--8<--
	disruption_py/settings/shot_ids_request.py:get_shot_ids_request_dict
	--8<--
	```
- A file path as a string with its suffix mapped to a ShotIdsRequest type in the `_file_suffix_to_shot_ids_request` dictionary:
	```python
	--8<--
	disruption_py/settings/shot_ids_request.py:file_suffix_to_shot_ids_request_dict
	--8<--
	```
- A dictionary mapping tokamak type strings to the desired shot ids request option for that tokamak.  E.g. `{'cmod': 'paper'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"

The following documents the support for shot ids requests:

::: disruption_py.settings.shot_ids_request
    handler: python
	options:
	  heading: 2
	  members:
	  - ShotIdsRequestParams
	  - ShotIdsRequest

## Built-in Implemenations

::: disruption_py.settings.shot_ids_request.IncludedShotIdsRequest

::: disruption_py.settings.shot_ids_request.FileShotIdsRequest

::: disruption_py.settings.shot_ids_request.DatabaseShotIdsRequest