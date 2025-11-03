#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "${BASH_SOURCE[0]}")/setup.sh"

# virtual environment
# shellcheck source=/dev/null
OLD=$PS1
source "$DISPY_PYVERS_DIR/bin/activate"
export PS1="(DisPy:$DISPY_BRANCH-py$DISPY_PYVERS) $OLD"

