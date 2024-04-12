#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "${BASH_SOURCE[0]}")/setup.sh"

# virtual environment
# shellcheck source=/dev/null
source "$DISPY_PYVERS_DIR/bin/activate"
export PS1="(DisPy:${PS1#(}"

