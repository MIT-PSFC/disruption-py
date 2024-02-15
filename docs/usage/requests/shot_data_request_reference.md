

## Built-in implementations { .doc .doc-heading }
The following file defines the list of built-in shot data requests. Built-in parameter methods are members of these objects. To view the methods in these files, please see the [GitHub repository](https://github.com/MIT-PSFC/disruption-py).
```python
--8<--
disruption_py/shots/parameter_methods/built_in.py
--8<--
```

## Custom Implementations { .doc .doc-heading }

DisruptionPy allows users to add custom data retrieval methods to DisruptionPy by creating decorated methods inside of subclasses of `ShotDataRequest` and passing an instance of the subclass in the `shot_data_request` parameter of [`ShotSettings`][disruption_py.settings.ShotSettings]. Please see the documentation on [parameter methods][custom-parameter-methods] for more information about creating your own parameter methods.

::: disruption_py.settings.shot_data_request
    handler: python
	options:
	  heading_level: 3
	  show_root_heading: false
	  show_root_toc_entry: false
	  members:
	  - ShotDataRequest
	  - ShotDataRequestParams
