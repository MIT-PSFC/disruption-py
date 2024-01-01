A module for handling set times requests passed in the [`ShotSettings`][disruption_py.settings.ShotSettings] class. 
Set times requests are used by disruption_py to set the timebase for data retrieval from MDSPlus and any sql tables.

This module defines the abstract class [`SetTimesRequest`][disruption_py.settings.set_times_request.SetTimesRequest] that can have subclasses passed as the
`set_times_request` parameter to the `ShotSettings` class.
It also provides built_in classes and mappings to easily set the timebase for data retrievel for common use cases.

Currently, these are the options that can be passed to the `set_times_request` parameter in `ShotSettings`:

- An instance of a subclass of `SetTimesRequest`
- A string identifier in the `_set_times_request_mappings` dictionary:
```python
--8<--
disruption_py/settings/set_times_request.py:set_times_request_dict
--8<--
```
- A python list, numpy array, or pandas series (with the timebase as the values) that should be used as the times for the timebase. See [`ListSetTimesRequest`][disruption_py.settings.set_times_request.ListSetTimesRequest] for more details.
- A dictionary mapping tokamak type strings to the desired `SetTimesRequest` for that tokamak.  E.g. `{'cmod': 'efit'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"

The following documents the support for set times requests:

::: disruption_py.settings.set_times_request
    handler: python
	options:
	  heading_level: 2
	  members:
	  - SetTimesRequest
	  - SetTimesRequestParams

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.set_times_request
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^SetTimesRequest$"
		- "!^SetTimesRequestParams$"
