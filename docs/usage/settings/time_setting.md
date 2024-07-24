## Overview { .doc .doc-heading }
A module for handling time settings passed in the `RetrievalSettings` class. 
Set time settings used by DisruptionPy to set the timebase for data retrieval from MDSPlus and any SQL tables.

This module defines the abstract class [`TimeSetting`][disruption_py.settings.time_setting.TimeSetting] that can have subclasses passed as the
`time_setting` parameter to the `RetrievalSettings` class.
It also provides built-in classes and mappings to easily set the timebase for data retrievel for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed to the `time_setting` parameter in `RetrievalSettings`:

- An instance of a subclass of `TimeSetting`
- A string identifier in the `_time_setting_mappings` dictionary:

```python
--8<--
disruption_py/settings/time_setting.py:time_setting_dict
--8<--
```
- A Python list, NumPy array, or Pandas Series (with the timebase as the values) that should be used as the times for the timebase. See [`ListTimeSetting`][disruption_py.settings.time_setting.ListTimeSetting] for more details.
- A dictionary mapping tokamak type strings to the desired `TimeSetting` for that tokamak.  E.g. `{'cmod': 'efit'}`.

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.time_setting
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^TimeSetting"
		- "!^TimeSettingParams$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of time settings must inherit from the `TimeSetting` abstract class, implementing the abstract methods.

::: disruption_py.settings.time_setting
    handler: python
	options:
	  heading_level: 2
	  members:
	  - TimeSetting
	  - TimeSettingParams

