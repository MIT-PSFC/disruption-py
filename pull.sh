#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "${BASH_SOURCE[0]}")/setup.sh"

# folder
DISPY_LOG="$DISPY_DIR/logs/$(basename "${0%.sh}").$(date +%F)"
export DISPY_LOG=$DISPY_LOG
mkdir -p "$DISPY_LOG"

# poetry
which poetry
poetry self update \
> "$DISPY_LOG/poetry.log"

# for each repo folder
for FOLDER in "$DISPY_DIR"/repo/*
do

   # branch
   DISPY_BRANCH=$(basename "$FOLDER")
   export DISPY_BRANCH=$DISPY_BRANCH

   # log
   echo -e "\n$(date) :: $DISPY_BRANCH = $FOLDER"
   export LOG="$DISPY_LOG/$DISPY_BRANCH"
   mkdir -p "$LOG"
   pushd "$FOLDER" &> /dev/null || exit 10

   # reset
   {
      git status
      git reset --hard # "$DISPY_BRANCH"
      git clean -ffdx # -e .venv
      git pull
      git status
   } \
   2>&1 \
   > "$LOG/git.log"

   # check pyproject
   [[ -s pyproject.toml ]] || continue

   # for each python version
   for VENV in "$DISPY_DIR/venv/$DISPY_BRANCH-py"*
   do

      {

      # activate
      export DISPY_PYVERS="${VENV##*py}"
      source "$DISPY_DIR/repo/auto/activate.sh" || exit 11

      # log
      export LOG="$LOG/py$DISPY_PYVERS"
      mkdir -p "$LOG"
      echo -e "$(date) :: $DISPY_BRANCH @ py$DISPY_PYVERS"

      # env
      env \
      | grep -e ^DISPY_ -e PYTHONPATH \
      2>&1 \
      > "$LOG/env.log"

      # before
      poetry run pip list \
      2>&1 \
      > "$LOG/before.log"

      # upgrade
      pip install --upgrade pip \
      2>&1 \
      > "$LOG/upgrade.log"

      # install
      poetry install --with lab --with dev \
      2>&1 \
      > "$LOG/install.log"

      # after
      poetry run pip list \
      2>&1 \
      > "$LOG/after.log"

      # run
      poetry run python examples/efit.py \
      1> "$LOG/efit.out" \
      2> "$LOG/efit.err"

      # log
      echo -e "$(date) :: $DISPY_BRANCH @ py$DISPY_PYVERS = rc ${PIPESTATUS[0]}"

      # deactivate
      deactivate

      } &

   done

   wait
   popd &> /dev/null || exit 12

done \
2>&1 \
| tee "$DISPY_LOG/all.log"

