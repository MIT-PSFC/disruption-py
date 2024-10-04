
# DisruptionPy

An interoperable Python package for plasma disruption analysis and prediction using ML. 

## Background

A key element of plasma control systems (PCS) in tokamak reactors is the prediction and avoidance of disruptions, sudden losses of the thermal and magnetic energy stored within the plasma that can occur when tokamaks operate near regions of plasma instability or because of system malfunctions.
The energy released during disruptions can cause severe damage to plasma-facing components, limiting experimental operation or even the device lifetime.
This poses a serious challenge to next-step fusion experiments such as SPARC, which will have to operate near some of the limits of plasma stability to achieve its intended performance and will do so at for long and frequent intervals.
Previous work has shown the promise of machine-learning (ML) algorithms for disruption prediction in both DIII-D and EAST -- the Experimental Advanced Superconducting Tokamak in China -- PCS.
This is also due to the fact that fusion science currently lacks first-principle, theoretical solutions to fully predict and avoid disruptions. 

DisruptionPy is an open-source python package for training, updating, and evaluating algorithms for disruption prediction and avoidance that can be applied to Alcator C-Mod and DIII-D data, and can deploy models in DIII-D and EAST (TBD) PCSs.

## Overview

DisruptionPy makes it easy to retrieve tabular data from MDSplus databases efficiently.
Users can create their own methods and/or use built-in methods that retrieve and derive a variety of important parameters from experimental data for disruption analysis.
These methods are run across all provided sets of discharges (or shot ids), outputting tabular data in customizable formats.

## Project layout

```python
disruption_py/ # source code
docs/ # documentation
examples/ # example workflows
scripts/ # miscellaneous scripts
tests/ # automated testing
```

The original Matlab scripts are now stored in the `matlab` [protected branch](https://github.com/MIT-PSFC/disruption-py/tree/matlab).

## Installation

DisruptionPy is now open-source and [available at PyPI](https://pypi.org/project/disruption-py/)!
For standard installations, please follow the usual way:

```bash
# if you use poetry:
poetry add disruption-py

# if you use pip:
pip install disruption-py
```

For custom installations, please refer to our [Installation guide](docs/INSTALL.md).

## Getting Started

Please see the project [quickstart](https://mit-psfc.github.io/disruption-py/quickstart/usage_quickstart/).

## Issues

If you have an issue please crate an issue on the GitHub repository

## Development

Please create a pull request if you have something to contribute!
