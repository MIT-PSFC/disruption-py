## Overview { .doc .doc-heading }
A module for handling the cache setting passed in the [`ShotSettings`][disruption_py.settings.ShotSettings] class. 
Cache settings are used by disruption_py to get the data that has already been retrieved before data is retrieved 
from MDSplus. E.g. data already exists in the disruption_warnings sql table.

This module defines the abstract class [`CacheSetting`][disruption_py.settings.cache_setting.CacheSetting] that can have subclasses passed as the
`cache_setting` parameter to the `ShotSettings` class.
It also provides built_in classes and mappings to easily retrieve existing data for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed to the `cache_setting` parameter in `ShotSettings`:

- An instance of a subclass of `CacheSetting`
- A string identifier in the `_cache_setting_mappings` dictionary:
```python
--8<--
disruption_py/settings/cache_setting.py:cache_setting_dict
--8<--
```
- A dictionary mapping tokamak type strings to the desired `CacheSetting` for that tokamak.  E.g. `{'cmod': 'sql'}`.
	<!-- --8<-- "disruption_py/utils/mappings/tokamak.py:allowed_tokamak_types_snippet" -->

## Built-in Implemenations { .doc .doc-heading }
::: disruption_py.settings.cache_setting
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^CacheSetting$"
		- "!^CacheSettingParams$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of cache setting must inherit from the `CacheSetting` abstract class, implementing the abstract methods.

::: disruption_py.settings.cache_setting
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		members:
		- CacheSetting
		- CacheSettingParams