# Getting Started with MAST Open Data & DisruptionPy

Welcome to DisruptionPy! This guide will help you install the software and run your first data retrieval workflows.

## What is DisruptionPy?

DisruptionPy is an open-source python package for efficient retrieval of experimental fusion data. It enables fast database preparation for disruption analysis and machine learning applications across multiple tokamak devices including C-Mod, DIII-D, EAST, HBT-EP, and MAST.

In this guide, we will focus on using DisruptionPy with MAST Open Data, which does not require any special permissions to access.

---

## Installation

DisruptionPy is available on [PyPI](https://pypi.org/project/disruption-py/) and can be installed using modern Python package managers.

### Option 1: Installation with uv

[uv](https://docs.astral.sh/uv/) is a fast Python package manager. To install DisruptionPy:

```bash
uv tool install disruption-py
```

### Option 2: Installation with Poetry

[Poetry](https://python-poetry.org/) provides robust dependency management:

```bash
# Add disruption-py
poetry add disruption-py

# Activate the virtual environment
poetry run disruption-py
```

### Option 3: Installation with pip

For a simple pip installation:

```bash
# Install disruption-py
pipx install disruption-py
```

---

## Command Line Interface Examples

DisruptionPy provides a command-line interface for quick data retrieval without writing Python scripts.

### CLI Quick Start

Retrieve data for a MAST shot using the EFIT timebase:

```bash
disruption-py --tokamak mast --time-base efit --output my_data.csv 30421
```

**Options explained:**

- `--tokamak mast`: Specifies the MAST tokamak as the data source
- `--time-base efit`: Uses the EFIT reconstruction timebase for data interpolation
- `--output my_data.csv`: Saves retrieved data to a CSV file named `my_data.csv`
- `30421`: The shot number to retrieve data from


### CLI Example 1: EFIT Timebase with Specific Methods

Retrieve EFIT parameters on the EFIT timebase:

```bash
disruption-py \
  --tokamak mast \
  --time-base efit \
  --methods get_efit_parameters \
  --output mast_efit_data.csv \
  30421
```

### CLI Example 2: Disruption Warning Timebase

Use the disruption warning timebase:

```bash
disruption-py \
  --tokamak mast \
  --time-base disruption_warning \
  --output mast_disruption_warning.csv \
  30421
```

### CLI Example 3: Plasma Current (Ip) Timebase

Use the plasma current timebase:

```bash
disruption-py \
  --tokamak mast \
  --time-base ip \
  --output mast_ip_timebase.csv \
  30421
```

### CLI Example 4: Multiple Shots with Parallel Processing

Process multiple MAST shots in parallel:

```bash
disruption-py \
  --tokamak mast \
  --time-base efit \
  -m get_efit_parameters -m get_ip -m get_density \
  --processes 2 \
  --output mast_multiple_shots.csv \
  30420 30421 30422
```

 - The `-m`/`--methods` parameter specifies which methods to call on the backend object to retrieve data. Multiple `-m` flags can be used to extract different types of data (e.g., EFIT parameters, plasma current, and density) in a single run.

 - The `-p`/`--processes` option enables parallel processing to speed up data retrieval for multiple shots.

### CLI Example 5: Output to NetCDF Format

Save data in NetCDF format for xarray compatibility:

```bash
disruption-py \
  --tokamak mast \
  --time-base efit \
  --output mast_data.nc \
  30421
```

### CLI Help

For a full list of available options:

```bash
disruption-py --help
```

## Python Usage Examples

### Quick Start Example

Here's a minimal example to retrieve data from a tokamak shot:

```python
from disruption_py.workflow import get_shots_data
from disruption_py.settings import RetrievalSettings

# Define retrieval settings
retrieval_settings = RetrievalSettings(
    time_setting="efit",  # Use EFIT timebase
)

# Retrieve data for MAST shot
shot_data = get_shots_data(
    tokamak="mast",
    shotlist_setting=[30420],
    retrieval_settings=retrieval_settings,
    output_setting="dataset",
)

print(shot_data)
```

---

## Working with Different Time Bases

One of DisruptionPy's key features is the ability to retrieve data on different timebases. The timebase determines which time points are used for data interpolation.

### Available Time Settings

DisruptionPy provides several built-in time settings:

| Time Setting | Description |
|-------------|-------------|
| `"efit"` | Uses EFIT reconstruction timebase |
| `"disruption_warning"` | Uses disruption warning timebase (machine-specific) |
| `"ip"` | Uses plasma current (Ip) timebase |
| Custom list/array | Use your own time points |

### Example 1: Using EFIT Timebase

EFIT provides equilibrium reconstruction at specific time points:

```python
from disruption_py.workflow import get_shots_data
from disruption_py.settings import RetrievalSettings

# Use EFIT timebase
retrieval_settings = RetrievalSettings(
    time_setting="efit",
    run_methods=["get_efit_parameters"],  # Retrieve EFIT parameters
)

shot_data = get_shots_data(
    tokamak="mast",
    shotlist_setting=[30420],
    retrieval_settings=retrieval_settings,
    output_setting="dataset",
)

print(f"Retrieved {len(shot_data.time)} time points")
print(f"Time range: {shot_data.time.min():.3f} - {shot_data.time.max():.3f} s")
```

### Example 2: Using Disruption Warning Timebase

The disruption warning timebase is optimized for disruption prediction:

```python
from disruption_py.workflow import get_shots_data
from disruption_py.settings import RetrievalSettings

# Use disruption warning timebase
retrieval_settings = RetrievalSettings(
    time_setting="disruption_warning",
)

shot_data = get_shots_data(
    tokamak="mast",
    shotlist_setting=[30420],
    retrieval_settings=retrieval_settings,
)

print(shot_data)
```

### Example 3: Using a Custom Time Array

You can provide your own time points as a list or NumPy array:

```python
import numpy as np
from disruption_py.workflow import get_shots_data
from disruption_py.settings import RetrievalSettings

# Define custom time points (in seconds)
custom_times = np.linspace(0.1, 0.5, 100)  # 100 points from 0.1 to 0.5 seconds

retrieval_settings = RetrievalSettings(
    time_setting=custom_times,
)

shot_data = get_shots_data(
    tokamak="mast",
    shotlist_setting=[30420],
    retrieval_settings=retrieval_settings,
)

print(f"Retrieved data at {len(shot_data.time)} custom time points")
```

### Example 4: Using Plasma Current (Ip) Timebase

You can also use the plasma current timebase:

```python
from disruption_py.workflow import get_shots_data
from disruption_py.settings import RetrievalSettings

# Use plasma current timebase
retrieval_settings = RetrievalSettings(
    time_setting="ip",
)

# Process MAST shot
shot_data = get_shots_data(
    tokamak="mast",
    shotlist_setting=[30420],
    retrieval_settings=retrieval_settings,
)

print(f"Retrieved {len(shot_data.time)} time points from Ip timebase")
```

---

## Getting Help

- **Documentation**: Visit the [full documentation](https://mit-psfc.github.io/disruption-py/)
- **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/MIT-PSFC/disruption-py/issues)
- **Contact**: Reach out to the development team at disruption-py@lists.psfc.mit.edu

---
