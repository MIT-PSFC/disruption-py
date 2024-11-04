
# DisruptionPy

#### An open-source physics-based Scientific Framework for Disruption Analysis of Fusion Plasmas for AI/ML applications

[![Workflow: Lint](https://github.com/MIT-PSFC/disruption-py/actions/workflows/lint.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/lint.yml)
[![Workflow: Tests](https://github.com/MIT-PSFC/disruption-py/actions/workflows/tests.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/tests.yml)
[![Workflow: Build](https://github.com/MIT-PSFC/disruption-py/actions/workflows/build.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/build.yml)
[![Workflow: Docs](https://github.com/MIT-PSFC/disruption-py/actions/workflows/docs.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/docs.yml)
[![Workflow: Dependabot](https://img.shields.io/badge/Dependabot-enabled-34d058?logo=github)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/dependabot/dependabot-updates)
[![Workflow: Stale](https://img.shields.io/badge/Stale%20bot-enabled-34d058?logo=github)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/stale.yml)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)
[![Linting: ruff](https://img.shields.io/badge/linting-ruff-purple)](https://github.com/astral-sh/ruff)
[![Linting: shellcheck](https://img.shields.io/badge/linting-shellcheck-lightgreen)](https://github.com/koalaman/shellcheck)
[![Linting: yamllint](https://img.shields.io/badge/linting-yamllint-lightblue)](https://github.com/adrienverge/yamllint)
[![Testing: pytest](https://img.shields.io/badge/testing-pytest-red)](https://github.com/pylint-dev/pylint-pytest)

[![Supported versions](https://img.shields.io/pypi/pyversions/disruption-py)](pyproject.toml)
[![Stats: downloads](https://static.pepy.tech/badge/disruption-py)](https://pepy.tech/project/disruption-py)
[![Available: PyPI](https://img.shields.io/pypi/v/disruption-py.svg)](https://pypi.org/project/disruption-py/)
[![Available: Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.13935223.svg)](https://doi.org/10.5281/zenodo.13935223)
[![License: MIT](https://img.shields.io/pypi/l/disruption-py?color=750014)](LICENSE)

## Concept

DisruptionPy is an open-source Scientific Python package for fast retrieval of experimental Fusion data from [MDSplus](https://www.mdsplus.org/) servers.
The library allows an efficient database preparation for downstream analysis and/or ML model development for disruption studies.
At present, the main supported machines are [Alcator C-Mod](https://en.wikipedia.org/wiki/Alcator_C-Mod) and [DIII-D](https://en.wikipedia.org/wiki/DIII-D_(tokamak)).


## Overview
### Background
A key element to ensure steady state operations in magnetically confined tokamak devices is the prediction and avoidance of disruptions.
These are sudden losses of the thermal and magnetic energy stored within the plasma, which can occur when tokamaks operate near stability boundaries or because of hardware anomalies.
The energy stored in the plasma and released during disruptions over milliseconds can cause severe damage to plasma-facing components, limiting experimental operations and the device's lifespan [[1](https://www.tandfonline.com/doi/full/10.1080/15361055.2023.2229675)].
Disruptions still pose a serious challenge to next-generation fusion devices such as ITER or SPARC, which will have to operate near some of the limits of plasma stability to achieve intended performance and will do so at for long and frequent intervals.
Fusion science currently lacks first-principle, theoretical solutions to fully predict and avoid disruptions. 
However, previous work [[2](https://doi.org/10.1088/1741-4326/ab28bf),[3](https://doi.org/10.1088/1741-4326/abf74d)] has shown the usefulness of machine-learning (ML) algorithms for disruption prevention for both DIII-D and EAST -- the Experimental Advanced Superconducting Tokamak in China -- operations.
DisruptionPy provides a standardized analysis pipeline across different fusion devices to build ML-ready datasets.

### Workflow
DisruptionPy makes it easy to retrieve experimental data from MDSplus [[4](https://www.mdsplus.org/index.php/Introduction)] fusion repositories efficiently.
Users can create their own routines and/or use built-in ones that retrieve and derive a variety of important signals from experimental data for disruption analysis.
These routines are then interpolated on a requested timebase across the specified set of plasma discharges (or shots) to assemble a dataset and save it under a variety of available formats.

< insert workflow image here > TODO

### Acknowledgments
DisruptionPy was fully revamped under the "Open and FAIR Fusion for Machine Learning Applications" [[5](https://crea-psfc.github.io/open-fair-fusion/)], with work partly funded by DOE ... TODO

### References
TODO

## Repository layout

Notable branches:

- `main`, the [stable branch](https://github.com/MIT-PSFC/disruption-py/tree/main),
- `dev`, the [development branch](https://github.com/MIT-PSFC/disruption-py/tree/dev),
- `matlab`, the [historical branch](https://github.com/MIT-PSFC/disruption-py/tree/matlab).


## Project layout

Brief description of the folders in our project:

- `disruption_py/`, package source code,
- `docs/`, documentation sources,
- `drafts/`, experimental scripts,
- `examples/`, example workflows,
- `scripts/`, miscellaneous scripts,
- `tests/`, testing workflows.


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

Please see the [project quickstart](https://mit-psfc.github.io/disruption-py/quickstart/usage_quickstart/).


## Contributing

> [!IMPORTANT]
> Make sure you refer to the latest version of our [development branch](https://github.com/MIT-PSFC/disruption-py/tree/dev)!

- If you encounter any problems, please [create a new issue](https://github.com/MIT-PSFC/disruption-py/issues/new).
- If you would like to contribute, please [submit a pull request](https://github.com/MIT-PSFC/disruption-py/compare/dev...).
- If you have general questions, please [start a new discussion](https://github.com/MIT-PSFC/disruption-py/discussions/new?category=q-a).
