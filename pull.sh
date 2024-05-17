#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "${BASH_SOURCE[0]}")/setup.sh"

# folder
TODAY=$(date +%F)
DISPY_LOG="$DISPY_DIR/logs/$(basename "${0%.sh}").$TODAY"
export DISPY_LOG=$DISPY_LOG
mkdir -p "$DISPY_LOG"

# poetry
which poetry
poetry self update \
> "$DISPY_LOG/poetry.log"

# API token
GAPI="https://api.github.com/repos/mit-psfc/disruption-py"
AUTH="$(cat "/home/$USER/.gh_pat")"

# for each repo folder
for FOLDER in "$DISPY_DIR"/repo/*
do

   # branch
   DISPY_BRANCH=$(basename "$FOLDER")
   export DISPY_BRANCH=$DISPY_BRANCH

   # args
   if [[ $# -ge 1 ]] && [[ "$1" != "$DISPY_BRANCH" ]]
   then
      continue
   fi

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

   # read statuses
   SHA=
   if [[ $DISPY_BRANCH =~ ^main ]] || [[ $DISPY_BRANCH =~ ^dev ]]
   then
      SHA=$(git rev-parse HEAD)
      curl -s \
      -H "$AUTH" \
      -D "$LOG/sha.txt" \
      -o "$LOG/sha.json" \
      "$GAPI/commits/$SHA/statuses" \
      2>&1 \
      > "$LOG/curl.log"
   fi

   # for each python version
   for VENV in "$DISPY_DIR/venv/$DISPY_BRANCH-py"*
   do

      # args
      if [[ $# -ge 2 ]] && [[ "$2" != "${VENV##*py}" ]]
      then
         continue
      fi

      {

      TSTART=$SECONDS

      # read status
      STATUS="Install / $(basename "$VENV") @ ${HOSTNAME%-*}"
      STATE=
      if [[ -n "$SHA" ]] && [[ -s "$LOG/sha.json" ]]
      then
         STATE=$(jq -r "sort_by(.updated_at) | map(select(.context == \"$STATUS\")) | last.state" "$LOG/sha.json" 2> /dev/null)
         [[ "$STATE" = "null" ]] && STATE=
      fi

      # activate
      export DISPY_PYVERS="${VENV##*py}"
      source "$DISPY_DIR/repo/auto/activate.sh" || exit 11

      # log
      export LOG="$LOG/py$DISPY_PYVERS"
      mkdir -p "$LOG"
      echo -e "$(date) :: $DISPY_BRANCH @ py$DISPY_PYVERS${STATE:+ = }$STATE"

      # done
      [[ "$STATE" = "success" ]] && exit 0

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

      # fast test
      export GITHUB_ACTIONS=1
      poetry run pytest -v tests \
      1> "$LOG/test.out" \
      2> "$LOG/test.err"
      RC=${PIPESTATUS[0]}

      # return code
      [[ $RC -eq 0 ]] && STATE=success || STATE=failure

      # log
      echo -e "$(date) :: $DISPY_BRANCH @ py$DISPY_PYVERS = rc $RC"

      # deactivate
      deactivate

      # status
      [[ -z "$SHA" ]] && exit 0
      echo "{\"state\":\"$STATE\",\"description\":\"Updated on $TODAY in $((SECONDS-TSTART)) s\",\"context\":\"$STATUS\"}" \
      | tee "$LOG/data.json" \
      | curl -s \
         -X POST \
         -H "$AUTH" \
         -D "$LOG/write.txt" \
         -o "$LOG/write.json" \
         "$GAPI/commits/$SHA/statuses" \
         -d @- \
         2>&1 \
         >> "$LOG/curl.log"

      } &

      wait # DEBUG

   done

   wait
   popd &> /dev/null || exit 12

done \
2>&1 \
| tee "$DISPY_LOG/all.log"

