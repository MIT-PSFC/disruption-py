## Custom Parameter Methods { .doc .doc-heading }

Users of disruption_py can create their own custom parameter methods, by adding decorators to methods in a subclass of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest]. Included methods with the `parameter_cached_method` decorator will have there results output alongside the results from methods in the respective `Shot` classes. See [`parameter_cached_method`][disruption_py.utils.method_caching.parameter_cached_method] for more details.

Decorated methods will be called with the following arguments:
TODO

::: disruption_py.settings.shot_data_request
    handler: python
	options:
	  heading_level: 3
	  show_root_heading: false
	  show_root_toc_entry: false
	  members:
	  - ShotDataRequest
	  - ShotDataRequestParams

## Decorators { .doc .doc-heading }

Methods inside of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest] subclasses can be decorated with the following 
decorators:

::: disruption_py.utils.method_caching
    handler: python
	options:
	  heading_level: 4
	  show_source: false
	  show_root_heading: false
	  show_root_toc_entry: false
	  show_root_members_full_path: true
	  members:
	  - parameter_cached_method
	  - cached_method

### Functions as decorator arguments { .doc .doc-heading }

::: examples.cached_method_params_function.cached_method_params_function
    handler: python
	options:
	  heading_level: 4
	  show_source: false
	  show_root_heading: false
	  show_root_toc_entry: false