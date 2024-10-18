
# DisruptionPy

An open-source physics-based Scientific Framework for Disruption Analysis of Fusion Plasmas. 

## Background

A key element of plasma control systems (PCS) in tokamak reactors is the prediction and avoidance of disruptions, sudden losses of the thermal and magnetic energy stored within the plasma that can occur when tokamaks operate near regions of plasma instability or because of hardware malfunctions.
The energy stored in the plasma and released during disruptions can cause severe damage to plasma-facing components, limiting experimental operation and the device lifetime [1](https://www.tandfonline.com/doi/full/10.1080/15361055.2023.2229675).
This poses a serious challenge to next-step fusion experiments such as SPARC, which will have to operate near some of the limits of plasma stability to achieve intended performance and will do so at for long and frequent intervals.
Fusion science currently lacks first-principle, theoretical solutions to fully predict and avoid disruptions. 
However, previous work [2,3] has shown the usefulness of machine-learning (ML) algorithms for disruption prevention for both DIII-D and EAST -- the Experimental Advanced Superconducting Tokamak in China -- operations.

DisruptionPy is an open-source, interoperable python package for fast data retrieval of experimental data analysis from MDSplus repositories. 
The library also allows database development for downstream ML model training for disruption studies. 
Its current implementation is available for Alcator C-Mod and DIII-D data servers.

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

## References
[1] Maris, A. D., Wang, A., Rea, C., Granetz, R., & Marmar, E. (2023). The Impact of Disruptions on the Economics of a Tokamak Power Plant. Fusion Science and Technology, 80(5), 636–652. https://doi.org/10.1080/15361055.2023.2229675
[2] 
[3]
