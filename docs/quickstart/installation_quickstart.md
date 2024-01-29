# Installation  

## On Specfic Clusters

### Pre-requirements
In order to access the specific clusters, i.e. Alcator C-Mod or DIII-D, a user agreement must first be signed. A local host must be identified, Cristina Rea (<crea@psfc.mit.edu>) will assist with the logistics.

### Alcator C-Mod  

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
source dpy-venv/bin/activate
```

Next install the DisruptionPy package locally. Note that installation in developer mode with the `-e` flag is now unavailable because of the migration to poetry.
```bash
pip install ~/dpy-experimental/disruption-py # No --user needed if installing in a virtual env
```

Next install other necessary package contents using the helper script by running:
```bash
disruption_py_setup
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