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
- A file path as a string with its suffix mapped to a `OutputTypeRequest` type in the `_file_suffix_to_output_type_request` dictionary:
	```python
	--8<--
	disruption_py/settings/output_type_request.py:file_suffix_to_output_type_request_dict
	--8<--
	```
- A dictionary mapping tokamak type strings to the desired `OutputTypeRequest` for that tokamak.  E.g. `{'cmod': 'efit'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"
- A python list of any other output type request option that can be passed as the `output_type_request` parameter in `ShotSettings` (all options listed previously). See [`OutputTypeRequestList`][disruption_py.settings.output_type_request.OutputTypeRequestList] for more details.

The following documents the support for output type requests:

::: disruption_py.settings.output_type_request
    handler: python
	options:
	  heading_level: 2
	  members:
	  - OutputTypeRequest
	  - ResultOutputTypeRequestParams
	  - FinishOutputTypeRequestParams

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.output_type_request
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^OutputTypeRequest$"
		- "!^ResultOutputTypeRequestParams$"
		- "!^FinishOutputTypeRequestParams$"
