# Generating Datasets
generate_datasets.py provides disruption_py users the ability to easily generate dataset files for any subset of implemented shot parameters. Parameters can either be grabbed from the relevant SQL database, calculated directly from the original MDSPlus data, or use both where the user can define the default and fallback behavior. 
* NOTE: Currently, disruption_py can only guarantee that the timebase for each shot as well as parameters grabbed directly from MDSPlus will be correct. Any parameters calculated based on MDSPlus data may contain notable discrepancies from the SQL database.
::: scripts.generate_datasets
## Example Usage
### Generate a dataset using SQL table 
```
python3 generate_datasets.py --data_source=1 --shotlist=./data/example_shotlist.txt --efit_tree=efit01
```

### Generate a dataset directly from MDSPlus data
```
python3 generate_datasets.py --data_source=3 --shotlist=./data/example_shotlist.txt --efit_tree=efit01
```