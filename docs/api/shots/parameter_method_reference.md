## What is a Parameter Method? { .doc .doc-heading }
In DisruptionPy, Parameter methods are methods that produce tabular data in a standardized time base. Parameter methods must be instance, class, or static methods of a subclass of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest] and take a single argument params that is an instance of [`ShotDataRequestParams`][disruption_py.settings.shot_data_request.ShotDataRequest].

## Built-in Parameter Methods { .doc .doc-heading }
Built-in parameter methods are defined inside of the `disruption_py.shots.parameter_methods` package. All built in methods are included through the built-in shot data requests in the `disruption_py.shots.parameter_methods.built_in.DEFAULT_SHOT_DATA_REQUESTS` list (shown below), where the built-in parameter methods are members of the listed shot data requests. To view the specific methods that are listed please see the [repository](https://github.com/MIT-PSFC/disruption-py) on GitHub.
```python
--8<--
disruption_py/shots/parameter_methods/built_in.py
--8<--
```

## Custom Parameter Methods { .doc .doc-heading }
Users of disruption_py can create their own custom parameter methods by adding decorators to methods in a subclass of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest]. Instances of these classes can then be passed as the `shot_data_request` parameter in the [`ShotSettings`][disruption_py.settings.ShotSettings], and there results will be included alongside those returned by the built-in methods. See [`parameter_cached_method`][disruption_py.shots.helpers.method_caching.parameter_cached_method] for more details on decorators.

The steps for creating a custom parameter method are as follows:

1. Create a subclass of [`ShotDataRequest`][disruption_py.settings.shot_data_request.ShotDataRequest]
```python
from disruption_py.settings.shot_data_request import ShotDataRequest

class MyShotDataRequest(ShotDataRequest):
	...
```

2. Add an instance or class method to the subclass that takes an argument named `params` of type [`ShotDataRequestParams`][disruption_py.settings.shot_data_request.ShotDataRequestParams] and returns a pandas DataFrame. The method must be decorated with the `parameter_cached_method` decorator. The arguments passed to the decorator are important for disruption_py to run efficiently. See [`parameter_cached_method`][disruption_py.shots.helpers.method_caching.parameter_cached_method] for more details about available parameters.
```python
from disruption_py.settings.shot_data_request import ShotDataRequest
from disruption_py.shots.helpers.method_caching import parameter_cached_method

class MyShotDataRequest(ShotDataRequest):

	@parameter_cached_method(...)
	def ***_method(self, params: ShotDataRequestParams) -> pd.Dataframe:
		...
```

3. To retrieve data from MDSplus use the `shot_props` attribute of the `params` object. The `ShotProps` class has a number of useful attributes with the most useful being listed below. See the `ShotProps` class for more details.
    - `params.shot_props.shot_id`: the shot id of the shot for which data is being retrieved.
	- `params.shot_props.tree_manager`: a reference to the MDSplus tree manager for the shot. This object should be used to open MDSplus trees instead of the regular MDSplus `Tree` class as it ensures that trees are both efficiently reused and closed when they are no longer needed. See ['TreeManager'] for details on how to use the tree manager.
    - `params.shot_props.times`: the timebase of the shot for which data is being retrieved as a numpy array of times. A common development pattern is using the `params.shot.interpolation_method` method (defaults to `interp1`) to interpolate the retrieved/computed values to the desired timebase.
??? example "Shot Data Request Examples"

    === "Kappa Area Parameter in C-Mod"

        ```python
        --8<--
        examples/shot_data_request_docs.py:kappa_area_request_example
        --8<--
        ```

    === "Another example"

        No more examples yet, please share your shot_data_requests if you are willing!

## Running Parameter Methods { .doc .doc-heading }
Users can use a number of built-in parameter methods and/or create their own methods.

For a parameter method to be run after calling [`get_shots_data`][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] it must meet the following conditions:

1. The method must be a member of an instance of a subclass of `ShotDatRequest` that is either:
	- included inside of the `disruption_py.shots.parameter_methods.built_in.DEFAULT_SHOT_DATA_REQUESTS` list (built-in method)
	- included inside of the `shot_data_request` parameter of the shot settings.

	??? note "Including multiple instances of the same class"
		It is possible for multiple instances of the same class with different intitialization arguments to the `shot_data_request` parameter of the shot settings. This can allow for parameter methods to have a wide range in functionality. To prevent all methods with the same name from being run in a request, please see [using functions as decorator arguments][functions-as-decorator-arguments]

2. The method must have the `parameter_cached_method` decorator with its `tokamak` parameter either not set or set to the tokamak that you are retrieving data from.

	??? example "Example decorator"
		```python
		from disruption_py.settings.shot_data_request import ShotDataRequest
		from disruption_py.shots.helpers.method_caching import parameter_cached_method
		from disruption_py.utils.mappings.tokamak import Tokamak

		class MyShotDataRequest(ShotDataRequest):

			@parameter_cached_method(tokamak=Tokamak.CMOD, ...)
			def cmod_retrieve_data(self, params: ShotDataRequestParams) -> pd.Dataframe:
				...
		```

3. The method is included to run via either the `run_methods`, `run_tags`, or `run_columns` parameters of the shot settings.
    - To be included via `run_methods`, the method name must be listed inside of `run_methods`
	- To be included via `run_tags`, the method must have a tag listed in the `tags` parameter of the `parameter_cached_method` decorator that is inlcuded in `run_tags`
	- To be inlcuded via `run_columns`, the method must have a column list in the `columns` parameter of the `parameter_cached_method` decorator that is included in `run_columns`
