#!/bin/bash

# parameters
export DISPY_PYVERS=${DISPY_PYVERS:=3.10}
export DISPY_BRANCH=${DISPY_BRANCH:=main}

# reset
if [[ -d /usr/local/mfe/disruptions ]]
then

   # C-MOD workstations
   export PATH=/usr/local/bin:/usr/sbin:/usr/bin:/usr/local/cmod/bin:/opt/thinlinc/bin
   export MDSPLUS_DIR=/usr/local/mdsplus
   export DISPY_DIR=/usr/local/mfe/disruptions/disruption-py

elif [[ -d /fusion/projects/disruption_warning ]]
then

   # DIII-D cluster
   export PATH=/usr/bin:/usr/sbin
   export MDSPLUS_DIR=/fusion/usc/c8/opt/mdsplus/alpha/7.139.59
   export DISPY_DIR=/fusion/projects/disruption_warning/disruption-py

else

   exit 1

fi

unset PYTHONPATH

# mdsplus
export PATH=$PATH:$MDSPLUS_DIR/bin
export PYTHONPATH=$MDSPLUS_DIR/python:$PYTHONPATH
export LD_LIBRARY_PATH=$MDSPLUS_DIR/lib:$LD_LIBRARY_PATH

# poetry
export PATH=$DISPY_DIR/poetry/bin:$PATH

# disruption-py
export PATH=$DISPY_DIR/poetry/bin:$PATH
export PATH=$DISPY_DIR/repo/auto/bin:$PATH
export DISPY_BRANCH_DIR=$DISPY_DIR/repo/$DISPY_BRANCH
export DISPY_PYVERS_DIR=$DISPY_DIR/venv/$DISPY_BRANCH-py$DISPY_PYVERS
export PYTHONPATH=$DISPY_BRANCH_DIR/disruption_py:$PYTHONPATH

