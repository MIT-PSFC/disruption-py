# Generating Datasets
generate_datasets.py provides disruption_py users the ability to easily generate dataset files for any subset of implemented shot parameters. Parameters can either be grabbed from the relevant SQL database, calculated directly from the original MDSPlus data, or use both where the user can define the default and fallback behavior. 
* NOTE: Currently, disruption_py can only guarantee that the timebase for each shot as well as parameters grabbed directly from MDSPlus will be correct. Any parameters calculated based on MDSPlus data may contain notable discrepancies from the SQL database.
::: scripts.generate_datasets
    handler: python
## Example Usage
### Generate a dataset using shots in SQL table 
```
python3 generate_datasets.py --data_source=1 --shotlist=./data/example_shotlist.txt --efit_tree=efit01
```

### Generate a dataset directly from MDSPlus data
```
python3 generate_datasets.py --data_source=3 --shotlist=./data/example_shotlist.txt --efit_tree=efit01
```

## Argument List
.. argparse-md:: scripts.generate_datasets

|Optional arguments|Description|
|------------------|-----------|
    |-h, --help            |show argparser help message and exit|
    |--shotlist SHOTLIST   |Path to file specifying shotlist|
    |--feature_cols FEATURE_COLS | Either a file or comma-separated list of desired feature columns
    |--output_dir OUTPUT_DIR| Path to generated data.|
    |--timebase_signal TIMEBASE_SIGNAL| Signal whose timebase will be used as the unifying timebase of the dataset. (Current choices are 'ip', 'disruption_timebase', and 'flattop')|
    |--efit_tree EFIT_TREE| Name of efit tree to use for each shot. If left as None, the script will use the get_efit_tree method in database.py.
    |--data_source {0,1,2,3}| 0: Default to SQL database then MDSPlus.\n1: Default to MDSPlus then SQL database.\n2: SQL database only.\n3: MDSPlus only.|
    |--unique_id UNIQUE_ID | Unique identifier for the dataset. Used to name the output files.
    |--log LOG | By default, generate_datasets will log to commandline but if this argument is true it will log to a file in the output directory|
    |--log_level {0,1,2,3,4,5}| **Notset**:0 **Debug**:1,**Info**:2,**Warning**:3,**Error**:4 **Critical**:5|
