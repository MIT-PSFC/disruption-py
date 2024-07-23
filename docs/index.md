# DisruptionPy
An interoperable Python package for plasma disruption analysis and prediction using ML. 

## Background
A key element of plasma control systems (PCS) in tokamak reactors is the prediction and avoidance of disruptions, sudden losses of the thermal and magnetic energy stored within the plasma that can occur when tokamaks operate near regions of plasma instability or because of system malfunctions. The energy released during disruptions can cause severe damage to plasma-facing components, limiting experimental operation or even the device's lifetime. This poses a serious challenge to next-step fusion experiments such as SPARC, which will have to operate near some of the limits of plasma stability to achieve its intended performance and will do so for long and frequent intervals. Previous work has shown the promise of machine learning (ML) algorithms for disruption prediction in both DIII-D and EAST -- the Experimental Advanced Superconducting Tokamak in China -- PCS. ML algorithms are also promising because fusion science currently lacks first-principle, theoretical solutions to fully predict and avoid disruptions. 

DisruptionPy is an open-source Python package for training, updating, and evaluating algorithms for disruption prediction and avoidance that can be applied to Alcator C-Mod and DIII-D data, and can deploy models in DIII-D and EAST (TBD) PCSs.

## Overview
DisruptionPy makes it easy to retrieve tabular data from MDSplus databases efficiently. Users can create their own methods and/or use built-in methods that retrieve and derive a variety of important parameters from experimental data for disruption analysis. These methods are run across all provided sets of discharges (or shot ids), outputting tabular data in customizable formats.

## Project layout
```python
disruption_py  # Source code
    ├── core  # Machine-agnostic data processing and retrieval infrastructure
    │   ├── physics_method
    │   └── utils
    ├── data  # Plain text shot lists
    ├── io  # Data sources
    ├── machine  # Physics methods for calculating parameters for multiple tokamaks
    │   ├── cmod
    │   └── d3d
    └── settings  # Settings classes
docs # Mkdocs generated documentation
drafts # Experimental physics methods and other code that are a work-in-progress
examples # Example usage of DisruptionPy
tests # Automated testing for DisruptionPy
└── utils
```
The original matlab scripts can be found in the matlab branch of the [Github repo](https://github.com/MIT-PSFC/disruption-py/tree/matlab). 