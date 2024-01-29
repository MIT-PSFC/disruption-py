# Installation  

## On Specfic Clusters

### Pre-requirements
In order to access the specific clusters, i.e. Alcator C-Mod or DIII-D, a user agreement must first be signed. A local host must be identified, Cristina Rea (<crea@psfc.mit.edu>) will assist with the logistics.

#### Create a virtual environment
Optionally create a new folder and virtual environment and activate it
```bash
mkdir ~/dpy-projects
cd ~/dpy-projects
python -m venv ~/dpy-projects/dpy-venv
source ~/dpy-projects/dpy-venv/bin/activate
```

### Standard installation  
If you are unable to access the project on GitHub, please see troubleshooting.

1. Install DisruptionPy:
```bash
pip install git+ssh://git@github.com/MIT-PSFC/disruption-py.git@develop#egg=disruption_py
```
2. Do other necessary package contents using the built-in helper script by running:
```bash
disruption_py_setup
```

### Editable installation
If you would like to modify the code and use it as an editable package, you may follow these instructions. Note that this assumes the use
of poetry as the dependency manager for your project. If you are unable to use poetry as a dependency manager for your project, alternate solutions
should be considered.

, you can do so by first cloning the package:
1. Clone DisruptionPy:
```bash
git clone git@github.com:MIT-PSFC/disruption-py.git
git switch develop
```

2. Install DisruptionPy
```bash
poetry add --editable /path/to/disruption_py
```

## Troubleshooting

### Local Installation
First copy the project directory to your home directory:
```bash
cd ~
mkdir ~/dpy-experimental
rm -rf ~/dpy-experimental/disruption-py # if you have already copied this previously
cp -R /home/joshlor/disruption-py ~/dpy-experimental/disruption-py
```

Optionally create a virtual environment and activate it
```bash
cd ~/dpy-experimental
python -m venv ~/dpy-experimental/dpy-venv
source ~/dpy-experimental/dpy-venv/bin/activate
```

Next install the DisruptionPy package locally. Note that installation in developer mode with the `-e` flag is now unavailable because of the migration to poetry.
```bash
pip install ~/dpy-experimental/disruption-py # No --user needed if installing in a virtual env
```



## Locally
TBD

note: if you are having trouble with installation please try to use mferws2.

run `pip install disruption_py`
run `disruption_py_setup`


# Development

## Troubleshooting

### Stuck on Poetry Install
In terminal run one of:
- `export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring`
- `pyenv shell system` and then `python3 -m keyring --disable`