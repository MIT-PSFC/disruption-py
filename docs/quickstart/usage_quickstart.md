
## Scripts (Recommended)
Creating a script gives you the full functionality of DisruptionPy. 

### Examples
Use [`basic_example_1.py`](https://github.com/MIT-PSFC/disruption-py/blob/main/examples/basic_example_1.py) or [`basic_example_2.py`](https://github.com/MIT-PSFC/disruption-py/blob/main/examples/basic_example_2.py) for a simple way to get started. Check out [`all_defaults_example.py`](https://github.com/MIT-PSFC/disruption-py/blob/main/examples/all_defaults_example.py) to see the all the default settings for retrieving data. 

### Creating a DisruptionPy script
1. **Create the shot data retrieval settings**
	The retrieval settings allow you to specify the settings for retrieving data for a single shot. You can specify details like the columns you want filled by the physics methods and the timebase domain you want used. See [`RetrievalSettings`](disruption_py.settings.retrieval_settings.RetrievalSettings) for more information and a complete list of settings options.

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
	The [`get_shots_data`](/usage/workflow_reference#disruption_py.workflow.get_shots_data) method is the main entry point for retrieving data from DisruptionPy as well as setting up connections to the SQL database and the MDSplus server. See [`get_shots_data`](/usage/workflow_reference#disruption_py.workflow.get_shots_data) for more details.
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
		without having to connect to the associated SQL database. See the [example](https://github.com/MIT-PSFC/disruption-py/blob/main/examples/no_database.py) or simply change the code above to pass
		the dummy database's default inititalizer as the database initializer to the handler class.
		note: this will result in data retrieval being incorrect for parameter methods that depend on 
		data retrieved from the SQL table eg. `time_until_disrupt`
		```python
		from disruption_py.io.sql import DummyDatabase
		from disruption_py.workflow import get_shots_data
			shot_data = get_shots_data(
				...
				DummyDatabase.initializer
				...
			)
		```