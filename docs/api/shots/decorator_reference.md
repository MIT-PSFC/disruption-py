
## Decorators { .doc .doc-heading }

Methods inside of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest] subclasses can be decorated with the following decorators:

::: disruption_py.shots.helpers.method_caching
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

## Functions as decorator arguments { .doc .doc-heading }

Functions may be passed as some of the parameters of the `parameter_cached_method` and `cached_method` decorators. 
These functions are called at runtime before the decorated method is run. They allow for the dynamic determination of
parameters such as `used_trees` that may not be defined before runtime.

??? example "Example"
	```python
	--8<--
	examples/decorator_functions_docs.py:decorator_functions_example
	--8<--
	```


### Decorators Function structure
::: examples.decorator_function_docs.cached_method_params_function
    handler: python
	options:
	  heading_level: 4
	  show_source: false
	  show_root_heading: false
	  show_root_toc_entry: false