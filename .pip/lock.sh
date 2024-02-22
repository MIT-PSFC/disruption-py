#!/bin/bash
set -e

PIP=$(dirname "$0")

command -v poetry
poetry check || poetry lock --no-update

{
   echo main
   grep 'tool.poetry.group.*dependencies' "$PIP/../pyproject.toml" \
   | cut -d . -f4
} \
| while read -r GROUP
do
   poetry export --only "$GROUP" > "$PIP/$GROUP.txt"
done \
| xargs -I {} echo -n "--with {} " \
| xargs poetry export > "$PIP/all.txt"

git status "$PIP"
