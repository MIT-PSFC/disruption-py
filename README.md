# DisruptionPy
An interoperable Python package for plasma disruption analysis and prediction using ML. 

## Background
A key element of plasma control systems (PCS) in tokamak reactors is the prediction and avoidance of disruptions, sudden losses of the thermal and magnetic energy stored within the plasma that can occur when tokamaks operate near regions of plasma instability or because of system malfunctions. The energy released during  disruptions can cause severe damage to plasma-facing components, limiting experimental operation or even the device lifetime. This poses a serious challenge to next-step fusion experiments such as SPARC, which will have to operate near some of the limits of plasma stability to achieve its intended performance and will do so at for long and frequent intervals. Previous work has shown the promise of machine-learning (ML) algorithms for disruption prediction in both DIII-D and EAST -- the Experimental Advanced Superconducting Tokamak in China -- PCS. This is also due to the fact that fusion science currently lacks first-principle, theoretical solutions to fully predict and avoid disruptions. 

DisruptionPy is an open-source python package for training, updating, and evaluating algorithms for disruption prediction and avoidance that can be applied to Alcator C-Mod and DIII-D data, and can deploy models in DIII-D and EAST (TBD) PCSs.

## Overview
DisruptionPy makes it easy to retrieve tabular data from MDSplus databases efficiently. Users can create their own methods and/or use built-in methods that retrieve and derive a variety of important parameters from experimental data for disruption analysis. These methods are run across all provided sets of discharges (or shot ids), outputting tabular data in customizable formats.

## Project layout
```python
disruption_py # Source code
docs # Mkdocs generated documentation
iris_requirements # requirements.txt for D3D iris cluster
examples # Example usage of DisruptionPy
scripts # Scripts for various disruption_py supported workflows
test # Automated testing for DisruptionPy
```

The original Matlab scripts are now stored in the `matlab` [protected branch](https://github.com/MIT-PSFC/disruption-py/tree/matlab).

## Installation

### Pre-requirements
In order to access the specific clusters, i.e. Alcator C-Mod or DIII-D, a user agreement must first be signed. A local host must be identified, Cristina Rea (<crea@psfc.mit.edu>) will assist with the logistics.

### Standard installation 

**For CMod, if you are on the mfe workstations, it is highly recommended to be on mferws02, mferws03, or mferws04.**

1. [Optional] Create a new project directory
	```bash
	mkdir ~/dpy-projects
	cd ~/dpy-projects
	```

2. [Optional, but recommended] Create a virtual environment

	Create a new folder and virtual environment and activate it
	```bash
	python -m venv ~/dpy-projects/dpy-venv
	source ~/dpy-projects/dpy-venv/bin/activate
	```

3. Install DisruptionPy

	Run:
	```bash
	pip install git+ssh://git@github.com/MIT-PSFC/disruption-py.git@develop#egg=disruption_py
	```
	If you are unable to access the project on GitHub, please see [troubleshooting](https://mit-psfc.github.io/disruption-py/installation/#troubleshooting).

4. Do other necessary setup tasks using the built-in helper script by running:
	```bash
	disruption_py setup
	```

For more information and alternate installation methods please see the [project site](https://mit-psfc.github.io/disruption-py/installation/).

## Getting Started
Please see the project [quickstart](https://mit-psfc.github.io/disruption-py/quickstart/usage_quickstart/).

## Issues
If you have an issue plase crate an issue on the GitHub repository

## Development
Please create a pull request if you have something to contribute!