

## Built-in implementations { .doc .doc-heading }
The following file defines the list of built-in shot data requests. To view the methods in these files, please see the GitHub repository.
--8<--
disruption_py/shots/parameter_methods/built_in.py
--8<--

## Custom Implementations { .doc .doc-heading }

DisruptionPy allows users to custom data retrieval methods by created decorated methods inside of subclasses of `shot_data_request`. Please see the documentation on [parameter methods][custom-parameter-methods] more information.

::: disruption_py.settings.shot_data_request
    handler: python
	options:
	  heading_level: 3
	  show_root_heading: false
	  show_root_toc_entry: false
	  members:
	  - ShotDataRequest
	  - ShotDataRequestParams
