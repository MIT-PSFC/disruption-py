#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "${BASH_SOURCE[0]}")/setup.sh"

# temporary files
TMPD="/tmp/$USER/disruption-py"
TMPF="$TMPD/cronfile"
mkdir -p "$TMPD" || exit 1
chmod go-rwx "$TMPD"

# extract
crontab -l \
> "$TMPF.old"

# build
envsubst \
< cronfile \
> "$TMPF.new"

# hashsums
sha256sum -b "$TMPF".{old,new}

# compare
diff --color=always "$TMPF".{old,new}

# install
# crontab "$TMPF.new"

