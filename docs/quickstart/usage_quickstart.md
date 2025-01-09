
## Scripts (Recommended)
Creating a script gives you the full functionality of DisruptionPy. 

### Examples
For a simple way to get started, check out [`simple.py`](https://github.com/MIT-PSFC/disruption-py/blob/main/examples/simple.py) or [`defaults.py`](https://github.com/MIT-PSFC/disruption-py/blob/main/examples/defaults.py) to see the all the default settings for retrieving data. 

### Creating a DisruptionPy script
1. **Create the shot data retrieval settings**
	The retrieval settings allow you to specify the settings for retrieving data for a single shot. You can specify details like the columns you want filled by the physics methods and the timebase domain you want used. See [`RetrievalSettings`][disruption_py.settings.retrieval_settings] for more information and a complete list of settings options.

	```python
	from disruption_py.settings.retrieval_settings import RetrievalSettings

	retrieval_settings = RetrievalSettings(
		# Use the efit timebase when returning data
		time_setting="efit",
		# Run all available methods
		run_tags=["all"],
	)
	```

2. **Call `get_shots_data`.** 
	[`workflow.py`][disruption_py.workflow] is the main entrypoint for retrieving data from DisruptionPy. The [`get_shots_data`][disruption_py.workflow.get_shots_data] method will be the primary method you use to get data, but [`workflow.py`][disruption_py.workflow] can also help you set up connections to the SQL database and the MDSplus server. See [`workflow.py`][disruption_py.workflow] for more details.
	```python
	from disruption_py.settings.retrieval_settings import RetrievalSettings
	from disruption_py.workflow import get_shots_data
	retrieval_settings = RetrievalSettings(
		...
	)
	shot_data = get_shots_data(
		tokamak="cmod",
		# Retrieve data for the desired shots
		shotlist_setting=[1150805012, 1150805013, 1150805014],
		# Use the created retrieval_settings
		retrieval_settings=retrieval_settings,
		# Automatically stream retrieved data to a csv file by passing in a file path ending in .csv
		output_setting="data.csv",
		# Use a single process to retrieve the data
		num_processes=1,
	)
	```

	??? question "I don't want to connect to the SQL database."
		DisruptionPy provides a dummy database class that allows users to retrieve data from MDSplus
		without having to connect to the associated SQL database. Pass
		the dummy database's default initializer as the database initializer to the handler class.
		Note: this will result in data retrieval being incorrect for parameter methods that depend on 
		data retrieved from the SQL table eg. `time_until_disrupt`
		```python
		--8<--
		docs/examples/no_database.py
		--8<--
		```