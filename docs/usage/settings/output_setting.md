## Overview { .doc .doc-heading }
A module for handling the output setting passed in the [`get_shots_data`][disruption_py.workflow.get_shots_data] 
method. The output setting is used to handle the output of data from DisruptionPy as it is retrieved. This may include collecting all the data from a request and returning it as a list or streaming outputted data to a file as it is retrieved.

This module defines the abstract class [`OutputSetting`][disruption_py.settings.output_setting.OutputSetting] that can have subclasses passed as the
`output_setting` argument to the `get_shots_data` method.
It also provides built-in classes and mappings to easily set the output type for common use cases.

### Usage { .doc .doc-heading }
Currently, these are the options that can be passed as the `output_setting` argument to `get_shots_data`:

- An instance of a subclass of `OutputSetting`
- A string identifier in the `_output_setting_mappings` dictionary:
```python
--8<--
disruption_py/settings/output_setting.py:output_setting_dict
--8<--
```
- A file path as a string with its suffix mapped to an `OutputSetting` type in the `_file_suffix_to_output_setting` dictionary:
```python
--8<--
disruption_py/settings/output_setting.py:file_suffix_to_output_setting_dict
--8<--
```
- A dictionary mapping tokamak type strings to the desired `OutputSetting` for that tokamak.  E.g. `{'cmod': 'list'}`.

- A Python list of any other output type request option that can be passed as the `OutputSetting` argument to `get_shots_data` (all options listed previously). See [`OutputSettingList`][disruption_py.settings.output_setting.OutputSettingList] for more details.

## Built-in Implemenations { .doc .doc-heading }
::: disruption_py.settings.output_setting
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		filters:
		- "!^OutputSetting"
		- "!^CompleteOutputSettingParams$"

## Custom Implementations { .doc .doc-heading }
Custom implementations of output type requests must inherit from the `OutputTypeRequest` abstract class, implementing the abstract methods.

::: disruption_py.settings.output_setting
    handler: python
	options:
		show_root_heading: false
		show_root_toc_entry: false
		show_root_members_full_path: true
		members:
		- OutputSetting
		- OutputSettingParams
		- CompleteOutputSettingParams

