`cache_setting` parameter to the [`RetrievalSettings`][disruption_py.settings.retrieval_settings] class. Some subclasses and mappings covering common use cases for retrieving data are provided.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed to the `cache_setting` parameter in [`RetrievalSettings`][disruption_py.settings.retrieval_settings]:

- An instance of a subclass of [`CacheSetting`][disruption_py.settings.cache_setting.CacheSetting]
- A string identifier in the `_cache_setting_mappings` dictionary:
```python
--8<--
disruption_py/settings/cache_setting.py:cache_setting_dict
--8<--
```
- A dictionary mapping tokamak type strings to the desired [`CacheSetting`][disruption_py.settings.cache_setting.CacheSetting] for that tokamak.  E.g. `{'cmod': 'sql'}`.
--8<-- "disruption_py/machine/tokamak.py:allowed_tokamak_types_snippet"

## Built-in Implementations { .doc .doc-heading }
::: disruption_py.settings.cache_setting
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^CacheSetting$"
		- "!^CacheSettingParams$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of cache setting must inherit from the [`CacheSetting`][disruption_py.settings.cache_setting.CacheSetting] abstract class, implementing the abstract methods.

::: disruption_py.settings.cache_setting
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		members:
		- CacheSetting
		- CacheSettingParams