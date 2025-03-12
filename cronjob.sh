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

# build global
envsubst \
< cronfile \
> "$TMPF.new"

# build local
MYCRON="cronfile-$(echo "$HOSTNAME" | grep -o '^[a-z]*')"
if [[ -s "$MYCRON" ]]
then
   envsubst \
   < "$MYCRON" \
   >> "$TMPF.new"
fi

# hashsums
sha256sum -b "$TMPF".{old,new}

# compare
diff --color=always "$TMPF".{old,new}
diff -q "$TMPF".{old,new} && exit 0

# install
echo -en "\ninstall? y/[n] "
read -r ANSWER
if [[ "$ANSWER" = "y" ]]
then
   crontab "$TMPF.new"
fi
