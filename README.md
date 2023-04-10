# disruption-warning-db-workflow
!!! Work in progress !!!!  
# DisruptionPy

## Installation  
Installation is cluster-specific:
### CMOD:  
The following command will install the disruption warning package locally in develop mode. In develop mode, changes to the directory used for installation will be reflected in the installed package.   
```
pip3 install --user -e /home/hmturner/disruption_py # No --user needed if installing in a virtual env
```
### D3D:
First we need to load the proper modules:
```
module load python/3
module unload gcc-4.9.2
module load gcc7/default
```
(OPTIONAL) Create and activate a new virtual env. I named it disruptions but feel free to name it whatever you want. 
```
python3 -m venv disruptions
source disruptions/bin/activate
```
Next, because of the age of the iris cluster, we install a special list of dependency packages and their versions.
```
pip3 install -r iris_requirements.txt
```
Finally, we'll install the package. The following command will install the disruption warning package locally in develop mode. In develop mode, changes to the directory used for installation will be reflected in the installed package.   
```
pip3 install --user -e /fusion/projects/disruption_warning/disruption-warning-db-workflow/
```

NOTE(not true just yet): The directory for installation is not the one I use for daily development and I will only push changes to it that have been tested. 
## Getting Started
### D3D
[Example script](https://github.com/crea-psfc/disruption-warning-db-workflow/tree/d3d/scripts/example.py)
## Development
Feel free to leave a PR if you're using the library but need a new feature.
### D3D 
Matlab source folders:
1. /fusion/projects/disruption_warning/software/matlab_programs
2. /fusion/projects/disruption_warning/software/peaking_factors_d3d
3. /fusion/projects/disruption_warning/software/peaking_factors_d3d/recalc_bradial  
Repo location: /fusion/projects/disruption_warning/disruption-warning-db-workflow
### Tags
INFO: Need to ask an expert for a description or explanation 
TODO(optional priority) + comment: Coding change/fix that needs to happen at the given priority level

## Documentation
- [Design Document](https://probable-argument-b7b.notion.site/Workflow-Design-Document-a04529032bda4a999f42e75182a43258)
- [Thesis Proposal](https://www.overleaf.com/read/xyhqcgvzssqb)
- Roadmap(TODO) 
