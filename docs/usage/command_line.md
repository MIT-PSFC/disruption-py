
DisruptionPy offers a number of built-in scripts through its CLI to make the process of retrieving data from MDSplus easier.
To use the CLI, simply run `disruption_py **command**` from the command line (prepend `poetry run` to the command if you are inside of an environment managed by poetry).

## `disruption_py setup` { .doc .doc-heading }
::: disruption_py.cli.setup_script.setup_cmod
    handler: python
    options:
		separate_signature: false
		show_root_heading: false
		show_root_toc_entry: false
		show_signature: false
        filters: ["!^_[^_]"]

## `disruption_py evaluate` { .doc .doc-heading }

The `disruption_py evaluate` command evaluates the accuracy of parameter methods in DisruptionPy. It runs all parameter methods with the all tag that are included in `disruption_py` by comparing values retrieved from MDSplus using DisruptionPy to ground truth values.

When complete, the command prints a short report on the methods that have suceeded and failed. Success criteria is having more that 95% of results within 1% of the ground truth. 

!!! note

	For some tokamaks, ground truth values were generated using MatLab. As MATLAB, uses different numeric gradient methods than numpy there will always be numeric differences between generated data and ground truth data. For this command, the numpy method `np.gradient` is replaced with the equivelant method from MATLAB to provide more insightful results.

## `disruption_py run generate_datasets` { .doc .doc-heading }
The `disruption_py run generate_datasets` command provides a basic interface for using `disruption_py` from the command line.


### Example Usage { .doc .doc-heading }

#### Generate a dataset using shots in SQL table { .doc .doc-heading }
```bash
disruption_py run generate_datasets \
	--data_source=1 \
	--shotlist=./data/example_shotlist.txt \
	--efit_tree=efit01
```
#### Generate a dataset directly from MDSPlus data { .doc .doc-heading }
```bash
disruption_py run generate_datasets \
	--data_source=3 \
	--shotlist=./data/example_shotlist.txt \
	--efit_tree=efit01
```

!!! note

	You can also run the generate datasets command from the `scripts/generate_datasets.py` script if you have cloned DisruptionPy. For example,
	Run `python scripts/generate_datasets.py ***arguments**`

### Argument List
Last updated Feb 13, 2024

| Argument | Description |
|---|---|
| -h, --help | Show this help message and exit |
| --shotlist SHOTLIST | Path to file specifying shotlist. Default: None |
| --tokamak TOKAMAK | Tokamak to use for data source. Supports Alcator C-Mod ("cmod" for cmod). Default: None |
| --num_processes NUM | Number of processes for data retrieval. Default: 1 |
| --only_requested_columns  | Only create dataset with requested columns. Default: False |
| --feature_cols COLS | File or comma-separated list of feature columns. Similar to run columns in `ShotSettings`. Default: None |
| --output_dir DIR | Path to generated data. Default: './output/' |
| --timebase_signal SIGNAL | Signal for unifying timebase of dataset. Default: "efit" |
| --efit_tree TREE | efit tree for each shot. Uses `get_efit_tree` method if None. Default: "analysis" |
| --data_source {0,1} | Data source preference. 0: SQL then MDSPlus. 1: MDSPlus only. Default: 1 |
| --unique_id ID | Unique identifier for dataset. Default: generated_id() |
| --log LOG | Log to file in output directory if true. Default: False |
| --log_level {0-5} | Log level. 0: Notset, 1: Debug, 2: Info, 3: Warning, 4: Error, 5: Critical. Default: 2 |
| --label {binary,none} | Timestep disruption label. Supports only binary. Default: 'binary' |
| --run_methods METHODS | Methods to run. Default: None |
| --run_tags TAGS | Tags to run. Default: None |
| --filter FILTER | Run `filter_dataset` method on dataset. Default: True |
| --produce_train_test | Whether to produce train and test datasets. Default: False |

!!! note

	Please run `disruption_py run generate_datasets --help` for a complete and up to date set of arguments.
