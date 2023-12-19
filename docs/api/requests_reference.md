Requests provide an easy and customizable manner to interact with disruption_py.

# Existing Data Requests
A module for handling existing data requests within the disruption_py framework.

This module defines the abstract class ExistingDataRequest that can have subclasses passed as the
existing_data_request parameter to the [ShotSettings][disruption_py.settings.ShotSettings] class.
It also provides built_in classes and mappings to easily retrieve existing data for common use cases.

These are the options currently available:

- A subclass of `ExistingDataRequest` from the `disruption_py.settings.existing_data_request` module.
- A string identifier in the `_existing_data_request_mappings` dictionary, e.g. `'sql'`.
- A dictionary mapping tokamaks to the existing data request for that tokamak. 
	Dictionary values may be any other option listed here, e.g. `{'cmod': 'sql'}`.

::: disruption_py.settings.existing_data_request
    handler: python
	options:
      filters: []


# Shot Data Requests
::: disruption_py.settings.shot_data_requests
    handler: python