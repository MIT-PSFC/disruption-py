## Installation  

### Pre-requirements
In order to access the specific clusters, i.e. Alcator C-Mod or DIII-D, a user agreement must first be signed. A local host must be identified, Cristina Rea (<crea@psfc.mit.edu>) will assist with the logistics.

### Standard installation 

1. [Optional] Create a new project directory
	```bash
	mkdir ~/dpy-projects
	cd ~/dpy-projects
	```

2. [Optional, but recommended] Create a virtual environment

	Create a new folder and virtual environment and activate it
	```bash
	python -m venv ~/dpy-projects/dpy-venv
	source ~/dpy-projects/dpy-venv/bin/activate
	```

3. Install DisruptionPy

	Run:
	```bash
	pip install git+ssh://git@github.com/MIT-PSFC/disruption-py.git@develop#egg=disruption_py
	```
	If you are unable to access the project on GitHub, please see [troubleshooting][trouble-accessing-github].

4. Do other necessary setup tasks using the built-in helper script by running:
	```bash
	disruption_py setup
	```

### Editable installation

1. Check that you have poetry installed by running
```bash
poetry -V
```
If you do not have poetry installed, please follow the instructions [here](https://python-poetry.org/docs/#installation).

2. Clone DisruptionPy:
```bash
git clone git@github.com:MIT-PSFC/disruption-py.git
git switch develop
```

3. Install DisruptionPy
	
	You may follow any of the following options depending on the use case:

	1. Work within the DisruptionPy package and let poetry manage the environment (easier)

		Steps:

		1. Navigate to the root of the project directory of DisruptionPy
		2. Install the package and poetry will automatically setup a python environment:
			```bash
			poetry install
			```
		3. Run the setup script:
			```bash
			poetry run disruption_py setup
			```

		Now when using disruption_py prepend commands with `poetry run`. For instance, when running a script use `poetry run python **script.py**`, or when running the cli use `poetry run disruption_py **command**`.

	2. Use your own environment (gives more control)

		Steps:

		1. Activate your virtual environment
		2. Navigate to the root of the project directory of DisruptionPy
		3. Install DisruptionPy and other required dependencies by running:
			```bash
			poetry install
			```
		4. Run the setup script:
			```bash
			disruption_py setup
			```

		Now you can use the package as normal as long as your virtual environment is activated (you do not need to prepend commands with `poetry run`)

## Troubleshooting

### Issues with package versioning
If the machine that you are working on has an outdated version of python, you may be unable to install.

#### On CMod
On the mfe workstations use mferws02 or mferws03

### Trouble accessing GitHub
If you are unable to access the GitHub repository you can manually install the package. Note that you may be installing an older version of DisruptionPy.

#### On CMod
You can do this by running:
```bash
pip install /home/joshlor/disruption-py
```

### Stuck on Poetry Install
In terminal run one of:

- `export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring`
- `pyenv shell system` and then `python3 -m keyring --disable`