# Contributing
## Pre-requirements:
- User access to DIII-D or Alcator C-Mod servers 
- GitHub account
- Cristina Rea (crea@psfc.mit.edu) has added you to MIT-PSFC GitHub organization
## Setting up development environment
We recommend forking the repository and cloning it to your local folder on whatever cluster you're using. Install in edit mode following the instructions in the installation guide. You should then be all be set to begin development. You're also welcome to use the main repository but please make sure to create your own branch to work in.
## GitHub Issues
GitHub Issues is recommended to track bugs as well as feature requests and development. If you find a bug, please submit an issue on the repository. If you'd like a significant feature change/addition or are planning to implement one yourself, please submit an issue as well. This will allow for useful discussion and collaboration among contributors.
## Coding and Style Guidelines
In general, try to follow Google coding style guidelines for python. Every new method, class, module, or script should not be committed without a descriptive docstring.
### Docstring format

#### Module
#### Class
""" One-line description.

More in-depth description. This should describe the purpose of the class and its uses.



"""
#### Shot Method
""" One-line description.

More in-depth description. This should describe in more detail how parameters are calculated and the basic theoretical reasoning for the chosen method.

Parameters
----------

Returns
-------

Original Authors
----------------
- Yourself or whoever wrote the original routines you are converting

Sources
---------
- Any relevant scripts either relative paths if stored in this repo, or URLs/absolute paths if stored elsewhere
- Helpful papers or documents 
"""
### Useful Links
DIII-D Matlab source folders:    
- /fusion/projects/disruption_warning/software/matlab_programs  
- /fusion/projects/disruption_warning/software/peaking_factors_d3d  
- /fusion/projects/disruption_warning/software/peaking_factors_d3d/recalc_bradial    
Repo location: /fusion/projects/disruption_warning/disruption-warning-db-workflow
### Tags
QUESTION: Need to ask an expert for a description or explanation 
TODO(optional priority) + comment: Coding change/fix that needs to happen at the optionally-given given priority level

## Submitting Changes
Please submit all code changes as PRs, pull requests, from your forked repository or personal branch to either the main branch or the relevant cluster-specific branch: d3d, cmod or east. PRs will be reviewed within 1-3 days moving forward.  
## Easy ways to start contributing
- Check out TODOs and QUESTIONs in the codebase for potentially easy contributions
- Add missing documentation to functions and calculations you are knowledgeable about
- Add/fix documentaiton on this website
## Communication Channels(TODO: Link)
- Slack: https://mit-psfc.slack.com/archives/C051HKPHQM9
- Zoom meeting: https://mit.zoom.us/my/c.rea

## Important Notices and Warnings
While the project is currently private, **please be sure to include no confidential information such as usernames or passwords directly in the codebase**. The repository will eventually become public. Also, **please don't commit data from servers that you don't have express permission to move off their network.**
