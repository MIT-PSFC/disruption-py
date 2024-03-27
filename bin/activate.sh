#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "$0")/setup.sh"

# virtual environment
# shellcheck source=/dev/null
source "$DISPY_PYVERS_DIR/bin/activate"
export PS1="(DisPy:${PS1#(}"

