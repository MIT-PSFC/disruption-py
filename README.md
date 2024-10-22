
# DisruptionPy

#### An open-source physics-based Scientific Framework for Disruption Analysis of Fusion Plasmas for AI/ML applications

[![Workflow: Lint](https://github.com/MIT-PSFC/disruption-py/actions/workflows/lint.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/lint.yml)
[![Workflow: Tests](https://github.com/MIT-PSFC/disruption-py/actions/workflows/tests.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/tests.yml)
[![Workflow: Build](https://github.com/MIT-PSFC/disruption-py/actions/workflows/build.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/build.yml)
[![Workflow: Docs](https://github.com/MIT-PSFC/disruption-py/actions/workflows/docs.yml/badge.svg)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/docs.yml)
[![Workflow: Dependabot](https://img.shields.io/badge/Dependabot-enabled-33dd44?logo=github)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/dependabot/dependabot-updates)
[![Workflow: Stale](https://img.shields.io/badge/Stale%20bot-enabled-33dd44?logo=github)](https://github.com/MIT-PSFC/disruption-py/actions/workflows/stale.yml)

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

Please see the [project quickstart](https://mit-psfc.github.io/disruption-py/quickstart/usage_quickstart/).


## Contributing

- If you encounter any problems, please [create a new issue](https://github.com/MIT-PSFC/disruption-py/issues/new).
- If you would like to contribute, please [submit a pull request](https://github.com/MIT-PSFC/disruption-py/compare/dev...).
- If you have general questions, please [start a new discussion](https://github.com/MIT-PSFC/disruption-py/discussions/new?category=q-a).
