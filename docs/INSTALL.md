
# Installation

Two _public_ installations are currently maintained automatically:

- on the C-MOD MFE Workstations, and
- on the DIII-D Omega cluster.

A _private_ installation is possible on any machine, either on a personal laptop, or on any experimental cluster.

## Public installations

The MIT PSFC Disruption Studies Group hosts a `disruption-py` installation in a NFS folder.
The following steps apply to both C-MOD workstations and the DIII-D cluster.

Snippets for quick addition to a user's `.bashrc` file are provided in the subsections below.

#### Python

A recent Python 3.12 version is installed as a [miniconda](https://docs.anaconda.com/free/miniconda/) distribution.
It might be used directly from the `"$DISPY_DIR"/miniconda` folder.

#### Dependency management

We use [Poetry](https://python-poetry.org/) and [uv](https://docs.astral.sh/uv/) for dependency management.

The helper scripts automatically activate the appropriate virtual environment, so direct Poetry or uv usage is not needed.
They may be used directly from the `"$DISPY_DIR"/poetry` or `"$DISPY_DIR"/uv` folders.

#### Branches

Two branches are installed publicly and kept up to date:

- `main`, for stable workflows;
- `dev`, for fresh features.

The target branch can be controlled through the `DISPY_BRANCH` environment variable. 

#### Virtual environments

For each branch, a virtual environment based off Python 3.12 is available for usage.

The helper scripts will choose the appropriate virtual environment based on the `DISPY_BRANCH` environment variable.
It may also be used directly from the `"$DISPY_DIR"/venv` folder.

#### Setup and activation

A setup script will set all the required environment variables to ensure functionality and reproducibility.
It has to be _sourced_, rather than executed, through:

```bash
source "$DISPY_DIR"/repo/auto/setup.sh
```

More often, a user may choose to directly _activate_ the chosen virtual environment, through:

```bash
source "$DISPY_DIR"/repo/auto/activate.sh
# or the alias defined in our snippet
disruption-activate
```

The helper scripts rely on the user adopting the [Bash shell](https://www.gnu.org/software/bash/).

#### Execution

Even without setup or activation, execution can be invoked through `disruption-python`, which should work as a drop-in replacement for the usual `python` command.
Since `disruption-python` is an executable, its shorthand execution is made possible by the presence of the `"$DISPY_DIR"/repo/auto` folder within the `PATH` environment variable, as suggested in the installation-specific snippets below.

For example, one may begin a disruption-py interactive session through:

```bash
disruption-python -ic "import disruption_py"
```

Or execute a disruption-py-based script through:

```bash
disruption-python workflow.py
```

A script can therefore be executed seamlessly through `disruption-python` by specifying a custom shebang:

```
user@host:~$ head -n1 workflow.py
#!/usr/bin/env disruption-python

user@host:~$ chmod +x workflow.py

user@host:~$ ./workflow.py
```

All these helper scripts are subject to change -- we welcome any suggestion to make the process even smoother.

### C-MOD

Suggested snippet to be appended to the user's `~/.bashrc` file:

```bash
# disruption-py
export DISPY_DIR=/usr/local/mfe/disruptions/disruption-py
export DISPY_BRANCH=main # default. or dev  
export PATH=$PATH:$DISPY_DIR/repo/auto:$DISPY_DIR/poetry/bin:$DISPY_DIR/uv
alias disruption-activate='source "$DISPY_DIR"/repo/auto/activate.sh'
```

### DIII-D

Suggested snippet to be appended to the user's `~/.bashrc` file:

```bash
# disruption-py
export DISPY_DIR=/fusion/projects/disruption_warning/disruption-py
export DISPY_BRANCH=main # default. or dev  
export PATH=$PATH:$DISPY_DIR/repo/auto:$DISPY_DIR/poetry/bin:$DISPY_DIR/uv
alias disruption-activate='source "$DISPY_DIR"/repo/auto/activate.sh'
```

### EAST

Suggested snippet to be appended to the user's `~/.bashrc` file:

```bash
# disruption-py
export DISPY_DIR=/project/disruption-py
export DISPY_BRANCH=main # default. or dev
export PATH=$PATH:$DISPY_DIR/repo/auto:$DISPY_DIR/poetry/bin:$DISPY_DIR/uv
alias disruption-activate='source "$DISPY_DIR"/repo/auto/activate.sh'
```

## Private installation

As Free and Open-Source Software (FOSS), disruption-py can also be installed on any machine.
We currently provide an installation guide for Ubuntu-based boxes, but generic Unix machines or Windows systems should support similar or equivalent steps. 

### Pre-requisites

Disruption-py currently needs non-python software to be installed as a pre-requisite:

1. [MDSplus](https://www.mdsplus.org/): to connect to MDSplus data servers,
2. SQL drivers: to connect to SQL database servers.

MDSplus can be installed using their [installation guide](https://www.mdsplus.org/index.php/Downloads).

On Ubuntu-based systems, SQL drivers might be installed for example through the [Microsoft ODBC Driver](https://learn.microsoft.com/en-us/sql/connect/odbc/linux-mac/installing-the-microsoft-odbc-driver-for-sql-server?view=sql-server-ver16) `msodbcsql18` package, or [FreeTDS](https://www.freetds.org/) `{tds,unix}odbc` packages.

Note:

On C-MOD workstations, MDSplus is pre-installed system-wide but its configuration is only set up for the system python through a small path file (e.g.: `/usr/lib/python3/dist-packages/mdsplus.pth`).
For virtual environments, even those created off the system python, this small path file needs to be copied over in the `site-packages` folder (e.g.: `venv/lib/python3.10/site-packages/`) as it is _not inherited_ upon creation of the virtual environment.

Alternatively, one might adopt the more standard and more obvious solution of adding the equivalent path to the `PYTHONPATH` environment variable, which does get read by virtual environments and inherited by subshells.

### Requirements

Python is obviously a requirement for disruption-py -- please make sure you are running the desired Python version before creating a new installation, e.g.:

```bash
which python
```

Python requirements should be installed from the Poetry-native lockfile `poetry.lock` committed to the repository in the main folder.
If necessary, backward-compatible pip-style `requirements.txt` files can be produced through the [Poetry export command](https://python-poetry.org/docs/cli/#export).

### Virtual environments

We _strongly encourage_ users to create a specific virtual environment for disruption-py usage.
If using Poetry, one will be created automatically, e.g.:

```bash
poetry install --with dev
```
