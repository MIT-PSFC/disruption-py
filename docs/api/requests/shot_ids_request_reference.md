A module for handling shot ids requests passed in the [`get_shots_data``][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] 
function. Shot ids requests are used by disruption_py to get the shot ids of shots that should have data retrieved from MDSplus.

This module defines the abstract class [`ShotIdsRequest`][disruption_py.settings.shot_ids_request.ShotIdsRequest] that can have subclasses passed as the `shot_ids_request`
argument to the `get_shots_data` function.
It also provides built_in classes and mappings to easily define shot ids for common use cases.

Currently, these are the options that can be passed as the `shot_ids_request` argument to `get_shots_data`:

- An isntance of a subclass of `ShotIdsRequest`
- A single shot id as an `int` or `str`
- A python list of shot ids as any combination of `int` or `str`
- A dictionary key as a string from the built-in mappings to data files in the `_get_shot_ids_request_mappings` dictionary: 
	```python
	--8<--
	disruption_py/settings/shot_ids_request.py:get_shot_ids_request_dict
	--8<--
	```
- A file path as a string with its suffix mapped to a `ShotIdsRequest` type in the `_file_suffix_to_shot_ids_request` dictionary:
	```python
	--8<--
	disruption_py/settings/shot_ids_request.py:file_suffix_to_shot_ids_request_dict
	--8<--
	```
- A dictionary mapping tokamak type strings to the desired shot ids request option for that tokamak.  E.g. `{'cmod': 'paper'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"
- A python list of any other shot id request option that can be passed as the `shot_ids_request` argument to `get_shots_data` (all options listed previously). All designated shot numbers will be concatanated and any duplicates will be removed.

The following documents the support for shot ids requests:

::: disruption_py.settings.shot_ids_request
    handler: python
	options:
	  heading_level: 2
	  members:
	  - ShotIdsRequest
	  - ShotIdsRequestParams

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.shot_ids_request
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^ShotIdsRequestParams$"
		- "!^ShotIdsRequest$"