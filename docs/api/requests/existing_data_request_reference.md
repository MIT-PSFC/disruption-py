## Overview { .doc .doc-heading }
A module for handling existing data requests passed in the [`ShotSettings`][disruption_py.settings.ShotSettings] class. 
Existing data requests are used by disruption_py to get the data that has already been retrieved before data is retrieved 
from MDSplus. E.g. data already exists in the disruption_warnings sql table.

This module defines the abstract class [`ExistingDataRequest`][disruption_py.settings.existing_data_request.ExistingDataRequest] that can have subclasses passed as the
`existing_data_request` parameter to the `ShotSettings` class.
It also provides built_in classes and mappings to easily retrieve existing data for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed to the `existing_data_request` parameter in `ShotSettings`:

- An instance of a subclass of `ExistingDataRequest`
- A string identifier in the `_existing_data_request_mappings` dictionary:
```python
--8<--
disruption_py/settings/existing_data_request.py:existing_data_request_dict
--8<--
```
- A dictionary mapping tokamak type strings to the desired `ExistingDataRequest` for that tokamak.  E.g. `{'cmod': 'sql'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"

## Built-in Implemenations { .doc .doc-heading }
::: disruption_py.settings.existing_data_request
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^ExistingDataRequest$"
		- "!^ExistingDataRequestParams$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of existing data requests must inheiret from the `ExistingDataRequest` abstract class, implementing the abstract methods.

::: disruption_py.settings.existing_data_request
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		members:
		- ExistingDataRequest
		- ExistingDataRequestParams