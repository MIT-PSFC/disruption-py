#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "$0")/setup.sh"

# cwd
cd "$DISPY_DIR/runner" || exit 1

# run
./run.sh

