
## Scripts (Recommended)
Creating a scripts gives you the full functionality of disruption_py. 
### Creating a disruption-py script
1. **Create an instance of the desired handler class.** 
	The handler class manages the main entry points for data retrieval as well as setting up connections to the SQL database. See [`CModHandler`][disruption_py.handlers.cmod_handler.CModHandler] for more details.
	```python
	from disruption_py.handlers.cmod_handler import CModHandler

	cmod_handler = CModHandler()
	```

	??? question "I don't want to connect to the SQL database."
		DisruptionPy provides a dummy database class that allows users to retrieve data from MDSplus
		without having to connect to the associated SQL database. Simply change the code above to pass
		the dummy database's default inititalizer as the database initializer to the handler class.
		note: this will result in data retrieval being incorrect for parameter methods that depend on 
		data retrieved from the SQL table eg. `time_until_disrupt`
		```python
		from disruption_py.databases.dummy_database import DummyDatabase
		from disruption_py.handlers.cmod_handler import CModHandler

		cmod_handler = CModHandler(database_initializer=DummyDatabase.default)
		```

2. **Create the shot settings.** 
	The shot settings help to define the desired data retrieval methods and the desired output type for data retrieval for each shot. See [`ShotSettings`][disruption_py.settings.shot_settings.ShotSettings] for more details and a complete list of settings options.
	```python
	from disruption_py.handlers.cmod_handler import CModHandler
	from disruption_py.settings.shot_settings import ShotSettings

	cmod_handler = CModHandler()

	shot_settings = ShotSettings(
		# use the efit timebase when returning data 
		set_times_request="efit",
		
		# run all available methods
		run_tags=["all"],
	)
	```

3. **Call `get_shots_data`.** 
	The `get_shots_data` method is the main entry point for retrieving data from DisruptionPy. See [`get_shots_data`][disruption_py.handlers.cmod_handler.CModHandler.get_shots_data] for more details.
	```python
	from disruption_py.handlers.cmod_handler import CModHandler
	from disruption_py.settings.shot_settings import ShotSettings

	cmod_handler = CModHandler()
	shot_settings = ShotSettings(
		... 
	)

	shot_data = cmod_handler.get_shots_data(
		# retrieve data for the list of provided shot numbers
		shot_ids_request=[1150805012, 1150805013, 1150805014,],

		# use the created shot_settings
		shot_settings=shot_settings,

		# stream retrieved data to the csv file
		output_type_request="ip_data.csv", 

		# use a single process to retrieve the data
		num_processes = 1,
	)
	```

## Command Line

DisruptionPy offers a number of built-in scripts through its CLI to make the process of retrieving data from MDSplus easier.
To use the CLI, simply run `disruption_py **command**` from the command line (prepend `poetry run` to the command if you are inside of an environment managed by poetry).

The commands available are listed below. You may also use `disruption_py --help` for more details.

### generate_datasets
The standard command for generating a dataset using disruption_py, allow for the generation DPRF compatible datasets for training and inference. Currently only supports CMod data. Run `disruption_py run generate_datasets --help` for information on available arguments.
To use run:
```bash
disruption_py run generate_datasets
```
