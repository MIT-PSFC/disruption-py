A module for handling output type requests passed in the [`ShotSettings`][disruption_py.settings.ShotSettings] class. 
Output type requests are used to handle the output of data from disruption_py as it is retrieved. This may inlude collecting all the data from a request and returning it as a list or streaming outputted data to a file as it is retrieved.

This module defines the abstract class [`OutputTypeRequest`][disruption_py.settings.output_type_request.OutputTypeRequest] that can have subclasses passed as the
`output_type_request` parameter to the `ShotSettings` class.
It also provides built_in classes and mappings to easily set the output type for common use cases.

Currently, these are the options that can be passed to the `output_type_request` parameter in `ShotSettings`:

- An instance of a subclass of `OutputTypeRequest`
- A string identifier in the `_output_type_request_mappings` dictionary:
```python
--8<--
disruption_py/settings/output_type_request.py:output_type_request_dict
--8<--
```
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

::: disruption_py.settings.existing_data_request
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^SetTimesRequest$"
		- "!^SetTimesRequestParams$"
