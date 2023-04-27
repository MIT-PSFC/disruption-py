# Installation  
## Locally 
TBD
## On Specfic Clusters
### Pre-requirements:
In order to access the specific clusters, i.e. Alcator C-Mod or DIII-D, a user agreement must first be signed. A local host must be identified, Cristina Rea (<crea@psfc.mit.edu>) will assist with the logistics.

### Alcator C-Mod:  
The following command will install the DisruptionPy package locally in developer mode. In developer mode, changes to the directory used for installation will be reflected in the installed package.   
```
pip3 install --user -e /usr/local/mfe/disruptions/disruption_py # No --user needed if installing in a virtual env
```
### DIII-D:
DIII-D computational clusters are accessed via gateway server `cybele.gat.com`.
`iris` and more recently `saga` are DIII-D computational clusters. Directories are shared across the different clusters, but unix environment are different.
When working on `iris`, first we need to load the proper modules (ignore on `saga`):
```
module load python/3
module unload gcc-4.9.2
module load gcc7/default
```
(OPTIONAL) Create and activate a new virtual env. In the example below it is named "disruptions", but feel free to name it whatever you want. 
```
python3 -m venv disruptions
source disruptions/bin/activate
```
Next, because of the age of the `iris` cluster, we install a special list of dependency packages and their versions.
```
pip3 install -r iris_requirements.txt
```
Finally, we'll install the package. The following command will install the DisruptionPy package locally in developer mode. In developer mode, changes to the directory used for installation will be reflected in the installed package.   
```
pip3 install --user -e /fusion/projects/disruption_warning/disruption-warning-db-workflow/ # Don't use user tag if in virtualenv
```

NOTE: The directory used for installation is not the one dedicated to daily development. Developers will only push changes to it that have been tested. 
