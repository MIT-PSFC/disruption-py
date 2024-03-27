#!/bin/bash

export DISPY_BRANCH=test

# setup/activate
# shellcheck source=/dev/null
source "$(dirname "$0")/activate.sh"

# cwd
cd "$DISPY_DIR/runner" || exit 1

# run
./run.sh

