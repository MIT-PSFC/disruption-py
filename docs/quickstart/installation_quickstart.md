# Installation  

## On Specfic Clusters

### Pre-requirements
In order to access the specific clusters, i.e. Alcator C-Mod or DIII-D, a user agreement must first be signed. A local host must be identified, Cristina Rea (<crea@psfc.mit.edu>) will assist with the logistics.

### Alcator C-Mod  

The following commands will install the DisruptionPy package locally in developer mode. In developer mode, changes to the directory used for installation will be reflected in the installed package.

```
cd ~
mkdir ~/dpy-experimental
rm -rf ~/dpy-experimental/disruption-py # if you have already copied this previously
cp -R /home/joshlor/disruption-py ~/dpy-experimental/disruption-py
pip3 install --user -e ~/dpy-experimental/disruption-py # No --user needed if installing in a virtual env
```


## Locally
TBD

note: if you are having trouble with installation please try to use mferws2.

run `pip install disruption_py`
run `disruption_py setup`


# Development

## Troubleshooting

### Stuck on Poetry Install
In terminal run one of:
- `export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring`
- `pyenv shell system` and then `python3 -m keyring --disable`