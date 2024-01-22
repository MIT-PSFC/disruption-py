# Terminology

## Parameter Method

Parameter methods are methods that produce tabular data in the standardized time base.

# Method Running

## Method selection

For a parameter method to be run it must meet the following conditions:

1. It must be a method of a subclass of the `ShotDatRequest` class that is either included inside of `disruption_py.shots.parameter_methods.built_in.DEFAULT_SHOT_DATA_REQUESTS` or has an instance passed in the `shot_data_request` parameter of the shot settings.

2. It must have the `parameter_cached_method` decorator and the `tokamak` parameter either set to the running tokamak or not set.

3. The method is included to run via either the `run_methods`, `run_tags`, or `run_columns` parameters of the shot settings.
    - To be included via `run_methods`, the method name must be listed inside of `run_methods`
	- To be included via `run_tags`, the method must have a tag listed in the `tags` parameter of the `parameter_cached_method` decorator that is inlcuded in `run_tags`
	- To be inlcuded via `run_columns`, the method must have a column list in the `columns` parameter of the `parameter_cached_method` decorator that is included in `run_columns`

# Creating Parameter Methods

