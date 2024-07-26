## What is a Physics Method? { .doc .doc-heading }
In DisruptionPy, physics methods are methods that produce tabular data in a standardized time base. Physics methods must take a single argument `params` that is an instance of [`PhysicsMethodParams`][disruption_py.core.physics_method.params].

## Built-in Physics Methods { .doc .doc-heading }
While you can define your own, existing built-in physics methods are defined inside of the `disruption_py.machine` package.

For more information on available methods please see the built-in method documentation pages:

- [CMod Physics Methods](cmod_built_in_method_reference.md)
- [DIII-D Physics Methods](d3d_built_in_method_reference.md)

### Evaluating accuracy and reliability of built-in physics methods
DisruptionPy provides testing infrastructure to compare the new Python-based physics method calculations with their Matlab equivalent with `pytest tests/test_against_sql.py`. 

## Custom Physics Methods { .doc .doc-heading }
Users of DisruptionPy can create their own custom physics methods by adding the [`@physics_method`][disruption_py.core.physics_method.decorator.physics_method] decorator to a method. These custom physics methods can then be passed as the `custom_physics_methods` parameter in [`RetrievalSettings`][disruption_py.settings.retrieval_settings] and their results will be included alongside those returned by the built-in methods. See [Physics Method Decorators](decorator_reference.md) for more details on decorators.

### Physics methods structure

::: examples.custom_physics_method.decorated_physics_method
    handler: python
	options:
	  heading_level: 4
	  show_source: false
	  show_root_heading: false
	  show_root_toc_entry: false

### Walkthrough { .doc .doc-heading }
The steps for creating a custom physics method are as follows:

1. Create a function that takes an argument named `params` of type [`PhysicsMethodParams`][disruption_py.core.physics_method.params.PhysicsMethodParams] and returns a Python dictionary. The method must be decorated with the [`physics_method`][disruption_py.core.physics_method.decorator.physics_method] decorator. The arguments passed to the decorator are important for DisruptionPy to run efficiently. See [`physics_method`][disruption_py.core.physics_method.decorator.physics_method] for more details about available arguments.
```python
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.physics_method.decorator import physics_method

@physics_method(...)
def ***_method(params: PhysicsMethodParams) -> dict:
	...
```

2. To retrieve data from MDSplus use the `params` ([`PhysicsMethodParams`][disruption_py.core.physics_method.params.PhysicsMethodParams]) object. It contains many useful attributes, among which are the following:
    - `params.shot_id`: the shot id of the shot for which data is being retrieved.
	- `params.mds_conn`: a wrapper around the MDSplus connection for the process, that makes it easier to get data for a shot. See [`MDSConnection`][disruption_py.io.mds.MDSConnection] for details.
    - `params.times`: the timebase of the shot for which data is being retrieved as a NumPy array of times. A common development pattern is using the `params.interpolation_method` method (defaults to `interp1`) to interpolate the retrieved/computed values to the desired timebase.
??? example "Shot Data Request Examples"

    === "Kappa Area Parameter in C-Mod"

        ```python
        --8<--
		examples/custom_physics_method.py:kappa_area_request_example
        --8<--
        ```

!!! warning
	When two output columns have the same name, the column that appears in the final dataset is not guaranteed. This issue will be fixed in the near future.

## Running Physics Methods { .doc .doc-heading }
Users can use a number of built-in physics methods and/or create their own methods.

For a physics method to be run after calling [`get_shots_data`][disruption_py.workflow.get_shots_data] it must meet the following conditions:

1. The method must either:
	- be included inside of the `disruption_py.machine.method_holders.py` list (built-in method)
	- be included inside of the `custom_physics_methods` argument of [`RetrievalSettings`][disruption_py.settings.retrieval_settings] when getting shot data.

2. The method must have the `physics_method` decorator with its `tokamak` parameter either not set or set to the tokamak that you are retrieving data from.

3. The method is included to run via either the `run_methods`, `run_tags`, or `run_columns` parameters of the shot settings.
    - To be included via `run_methods`, the method name must be listed inside of `run_methods`
	- To be included via `run_tags`, the method must have a tag listed in the `tags` parameter of the `physics_method` decorator that is included in `run_tags`
	- To be included via `run_columns`, the method must have a column list in the `columns` parameter of the `physics_method` decorator that is included in `run_columns`


Once all designated methods have been collected, DisruptionPy optimizes their execution order to minimize resource usage by using the information supplied in the `physics_method` decorator. Once reordering is complete, the methods are run.


::: disruption_py.core.physics_method.params
	handler: python
	options:
	  heading_level: 2
	  members:
	   - PhysicsMethodParams
