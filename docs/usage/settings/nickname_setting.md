## Overview { .doc .doc-heading }
A module for handling the EFIT tree nickname passed in the `RetrievalSettings` class. 
DisruptionPy uses the nickname setting to determine which MDSplus EFIT tree to get data from.

This module defines the abstract class `NicknameSetting` that can have subclasses passed as the
`efit_nickname_setting` parameter to the `RetrievalSettings` class.
It also provides built-in classes and mappings to easily set the nickname for data retrieval for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed to the `efit_nickname_setting` parameter in `RetrievalSettings`:

- An instance of a subclass of `NicknameSetting`
- A string identifier in the `_nickname_setting_mappings` dictionary:
```python
--8<--
disruption_py/settings/nickname_setting.py:nickname_setting_keys
--8<--
```
- A dictionary mapping tokamak type strings to the desired `NicknameSetting` for that tokamak.  E.g. `{'cmod': 'efit'}`.
--8<-- "disruption_py/machine/tokamak.py:allowed_tokamak_types_snippet"

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.nickname_setting
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		members:
		- StaticNicknameSetting
        - DefaultNicknameSetting
        - DisruptionNicknameSetting

## Custom Implementations { .doc .doc-heading }
Custom implementations of a nickname setting must inherit from the `NicknameSetting` abstract class, implementing the abstract methods.

::: disruption_py.settings.nickname_setting
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
        members:
        - NicknameSetting
        - NicknameSettingParams
