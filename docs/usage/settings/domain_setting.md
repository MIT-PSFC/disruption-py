## Overview { .doc .doc-heading }
DisruptionPy uses the time domain setting to specify the portion of the shot (ie rampup, flattop, full shot) for which to retrieve data.

This module defines the abstract class [`DomainSetting`][disruption_py.settings.domain_setting.DomainSetting] that can have subclasses passed as the
`domain_setting` parameter to the [`RetrievalSettings`][disruption_py.settings.retrieval_settings] class. 
It also provides built-in classes and mappings to easily set the domain for data retrieval for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed to the `domain_setting` parameter in [`RetrievalSettings`][disruption_py.settings.retrieval_settings]:

- An instance of a subclass of [`DomainSetting`][disruption_py.settings.domain_setting.DomainSetting]
- A string identifier in the `_domain_setting_mappings` dictionary:
```python
--8<--
disruption_py/settings/domain_setting.py:domain_setting_dict
--8<--
```
- A dictionary mapping tokamak type strings to the desired [`DomainSetting`][disruption_py.settings.domain_setting.DomainSetting] for that tokamak.  E.g. `{'cmod': 'flattop'}`.
--8<-- "disruption_py/machine/tokamak.py:allowed_tokamak_types_snippet"

## Built-in Implemenations { .doc .doc-heading }

::: disruption_py.settings.domain_setting
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^DomainSetting"
		- "!^DomainSettingParams$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of a domain setting must inherit from the [`DomainSetting`][disruption_py.settings.domain_setting.DomainSetting] abstract class, implementing the abstract methods.

::: disruption_py.settings.domain_setting
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
        members:
        - DomainSetting
        - DomainSettingParams
