#!/bin/bash

# setup
# shellcheck source=/dev/null
source "$(dirname "${BASH_SOURCE[0]}")/setup.sh"

# folder
DISPY_LOG="$DISPY_DIR/logs/$(basename "${0%.sh}").$(date +%F.%s)"
export DISPY_LOG=$DISPY_LOG
mkdir -p "$DISPY_LOG"

# poetry
{
   which poetry
   poetry --version
   poetry self update
   poetry --version
} \
> "$DISPY_LOG/poetry.log" \
2>&1

# uv
{
   which uv
   uv --version
   uv self update
   uv --version
} \
> "$DISPY_LOG/uv.log" \
2>&1

# API token
GAPI="https://api.github.com/repos/mit-psfc/disruption-py"
AUTH="$(cat ~/.gh_pat)"

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
   elif [[ "$DISPY_BRANCH" == "test" ]]
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
      git fetch origin "$DISPY_BRANCH"
      git status
      git reset --hard "origin/$DISPY_BRANCH"
      git pull
      git clean -ffdx
      git status
   } \
   > "$LOG/git.log" \
   2>&1

   # check pyproject
   [[ -s pyproject.toml ]] || continue

   # read statuses
   SHA=
   JSON=
   if [[ $DISPY_BRANCH =~ ^main ]] || [[ $DISPY_BRANCH =~ ^dev ]] || [[ $DISPY_BRANCH =~ ^east ]]
   then
      SHA=$(git rev-parse HEAD)
      JSON="$LOG/read.json"
      curl -s \
      -H "$AUTH" \
      -D "$LOG/read.txt" \
      -o "$JSON" \
      "$GAPI/commits/$SHA/statuses" \
      > "$LOG/read.log" \
      2>&1
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

      # context
      if [[ "$GITHUB_ACTIONS" = "true" ]]
      then
         STATUS="Tests / pytest (${DISPY_TOKAMAK^^})"
      else
         STATUS="Deploy / $(echo "$HOSTNAME" | grep -o '^[a-z]*')"
      fi

      # read status
      STATE=
      if [[ -n "$SHA" ]] && [[ -s "$JSON" ]]
      then
         STATE=$(jq -r "sort_by(.updated_at) | map(select(.context == \"$STATUS\")) | last.state" "$JSON" 2> /dev/null)
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
      > "$LOG/env.log" \
      2>&1

      # before
      poetry run pip list \
      > "$LOG/before.log" \
      2>&1

      # upgrade
      pip install --upgrade pip \
      > "$LOG/upgrade.log" \
      2>&1

      # install
      poetry install --with lab --with dev \
      > "$LOG/install.log" \
      2>&1

      # after
      poetry run pip list \
      > "$LOG/after.log" \
      2>&1

      # full test
      make test \
      1> "$LOG/test.out" \
      2> "$LOG/test.err"
      RC=${PIPESTATUS[0]}

      # return code
      [[ $RC -eq 0 ]] && STATE=success || STATE=failure

      # log
      echo -e "$(date) :: $DISPY_BRANCH @ py$DISPY_PYVERS = rc $RC"

      # deactivate
      deactivate

      # clock
      TELAP=$((SECONDS-TSTART))
      if [[ $TELAP -lt 60 ]]
      then
         TMESS="${TELAP}s"
      else
         TMESS="$((TELAP/60))m$((TELAP%60))s"
      fi

      # status
      [[ -z "$SHA" ]] && exit 0
      echo "{\"state\":\"$STATE\",\"description\":\"Deployed on $(date +%F) in $TMESS\",\"context\":\"$STATUS\"}" \
      | tee "$LOG/data.json" \
      | curl -s \
         -X POST \
         -H "$AUTH" \
         -D "$LOG/write.txt" \
         -o "$LOG/write.json" \
         "$GAPI/commits/$SHA/statuses" \
         -d @- \
         > "$LOG/write.log" \
         2>&1

      } &

      wait # DEBUG

   done

   wait
   popd &> /dev/null || exit 12

done \
2>&1 \
| tee "$DISPY_LOG/all.log"

