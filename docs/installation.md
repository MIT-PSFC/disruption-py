# Installation  
## Locally 
TBD
## On Specfic Clusters
### CMOD:  
The following command will install the disruption warning package locally in develop mode. In develop mode, changes to the directory used for installation will be reflected in the installed package.   
```
pip3 install --user -e /usr/local/mfe/disruptions/disruption_py # No --user needed if installing in a virtual env
```
### D3D:
First we need to load the proper modules(ignore on saga cluster):
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
pip3 install --user -e /fusion/projects/disruption_warning/disruption-warning-db-workflow/ # Don't use user tag if in virtualenv
```

NOTE(not true just yet): The directory for installation is not the one I use for daily development and I will only push changes to it that have been tested. 