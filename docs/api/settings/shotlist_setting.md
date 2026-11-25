## Overview { .doc .doc-heading }
A module for handling shot ids passed in the [`get_shots_data`][disruption_py.workflow.get_shots_data] 
method. DisruptionPy will retrieve MDSplus data for those shot ids.

This module defines the abstract class [`ShotlistSetting`][disruption_py.settings.shotlist_setting.ShotlistSetting] that can have subclasses passed as the `shotlist_setting` argument to the `get_shots_data` method.
It also provides built-in classes and mappings to easily define shot ids for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed as the `shotlist_setting` argument to `get_shots_data`:

- An instance of a subclass of `ShotlistSetting`
- A single shot id as an `int` or `str`
- A Python list of shot ids as any combination of `int` or `str`
- A dictionary key as a string from the built-in mappings to data files in the `_get_shotlist_setting_mappings` dictionary: 
	```python
	--8<--
    disruption_py/settings/shotlist_setting.py:get_shotlist_setting_dict
	--8<--
	```
- A file path as a string with its suffix mapped to a `ShotlistSetting` type in the `_file_suffix_to_shotlist_setting` dictionary:
	```python
	--8<--
    disruption_py/settings/shotlist_setting.py:file_suffix_to_shotlist_setting_dict
	--8<--
	```
- A dictionary mapping tokamak type strings to the desired shot ids option for that tokamak.  E.g. `{'cmod': 'cmod_test'}`.
- A Python list of any other shot id request option that can be passed as the `shotlist_setting` argument to `get_shots_data` (all options listed previously). All designated shot numbers will be concatenated and any duplicates will be removed.

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.shotlist_setting
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^ShotIdsRequestParams$"
		- "!^ShotIdsRequest$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of shot id settings must inherit from the `ShotlistSetting` abstract class, implementing the abstract methods.

::: disruption_py.settings.shotlist_setting
    handler: python
	options:
	  heading_level: 2
	  members:
	  - ShotIdsRequest
	  - ShotIdsRequestParams
