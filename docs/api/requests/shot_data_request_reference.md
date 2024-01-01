## Custom Parameter Methods { .doc .doc-heading }

Users of disruption_py can create their own custom parameter methods, by adding decorators to methods in a subclass of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest]. Included methods with the `parameter_cached_method` decorator will have there results output alongside the results from methods in the respective `Shot` classes. See [`parameter_cached_method`][disruption_py.utils.method_caching.parameter_cached_method] for more details.

The steps for creating a custom parameter method are as follows:

1. Create a subclass of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest]
```python
from disruption_py.settings.shot_data_request import ShotDataRequest

class MyShotDataRequest(ShotDataRequest):
	...
```

2. Add an instance or class method to the subclass that takes an argument named `params` of type [`ShotDataRequestParams`][disruption_py.settings.shot_data_request.ShotDataRequestParams] and returns a pandas DataFrame. The method must be decorated with the `parameter_cached_method` decorator. The arguments passed to the decorator are important for disruption_py to run efficiently. See [`parameter_cached_method`][disruption_py.utils.method_caching.parameter_cached_method] for more details about available parameters.
```python
from disruption_py.settings.shot_data_request import ShotDataRequest
from disruption_py.utils.method_caching import parameter_cached_method

class MyShotDataRequest(ShotDataRequest):

	@parameter_cached_method(...)
	def ***_method(self, params: ShotDataRequestParams) -> pd.Dataframe:
		...
```

3. To retrieve data from MDSplus use the `shot` attribute of the `params` object. The shot object has a number of useful attributes with the most useful being listed below. See the [`Shot`][disruption_py.shots.shot.Shot] class for more details.
    - `params.shot.get_tree_manager()` returns a reference to the MDSplus tree manager for the shot. This object should be used to open MDSplus trees instead of the regular MDSplus `Tree` class as it ensures that trees are both efficiently reused and closed when they are no longer needed. See ['TreeManager'] for details on how to use the tree manager.
    - `params.shot.get_shot_id()` returns the shot id of the shot for which data is being retrieved.
    - `params.shot.get_times()` returns the timebase of the shot for which data is being retrieved as a numpy array of times. A common development pattern is using the `params.shot.interpolation_method` method (defaults to `interp1`) to interpolate the retrieved/computed values to the desired timebase.
??? example "Shot Data Request Examples"

    === "Kappa Area Parameter in C-Mod"

        ```python
        --8<--
        examples/shot_data_request_docs.py:kappa_area_request_example
        --8<--
        ```

    === "Another example"

        No more examples yet, please share your shot_data_requests if you are willing!

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

::: examples.shot_data_request_docs.cached_method_params_function
    handler: python
	options:
	  heading_level: 4
	  show_source: false
	  show_root_heading: false
	  show_root_toc_entry: false

## Class Reference { .doc .doc-heading }

::: disruption_py.settings.shot_data_request
    handler: python
	options:
	  heading_level: 3
	  show_root_heading: false
	  show_root_toc_entry: false
	  members:
	  - ShotDataRequest
	  - ShotDataRequestParams