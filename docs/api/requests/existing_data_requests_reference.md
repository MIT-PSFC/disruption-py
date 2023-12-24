# Existing Data Requests
A module for handling existing data requests passed in the [ShotSettings][disruption_py.settings.ShotSettings] class. 
Existing data requests are used by disruption_py to get the data that has already been retrieved before data is retrieved 
from MDSplus. E.g. data already exists in the disruption_warnings sql table.

This module defines the abstract class ExistingDataRequest that can have subclasses passed as the
existing_data_request parameter to the [ShotSettings][disruption_py.settings.ShotSettings] class.
It also provides built_in classes and mappings to easily retrieve existing data for common use cases.

Currently, these are the options that can be passed to the existing_data_request parameter in 
[ShotSettings][disruption_py.settings.ShotSettings]:

- A subclass of `ExistingDataRequest` from the `disruption_py.settings.existing_data_request` module.
- A string identifier in the `_existing_data_request_mappings` dictionary:
```python
--8<--
disruption_py/settings/existing_data_request.py:existing_data_request_dict
--8<--
```
- A dictionary mapping tokamak type strings to the desired `ExistingDataRequest` for that tokamak.  E.g. `{'cmod': 'sql'}`.
	--8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet"

::: disruption_py.settings.existing_data_request
    handler: python
	options:
      filters: []