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
      for N in {0..9}
      do
         git status
         git reset --hard
         git clean -ffdx
         git pull && break
         REV=$((1+N*2))
         echo "reset: HEAD~$REV"
         git reset --hard "HEAD~$REV"
      done
      git status
   } \
   2>&1 \
   > "$LOG/git.log"

   # check pyproject
   [[ -s pyproject.toml ]] || continue

   # read statuses
   SHA=
   if [[ $DISPY_BRANCH =~ ^main ]] || [[ $DISPY_BRANCH =~ ^dev ]] || [[ $DISPY_BRANCH =~ ^east ]]
   then
      SHA=$(git rev-parse HEAD)
      curl -s \
      -H "$AUTH" \
      -D "$LOG/read.txt" \
      -o "$LOG/read.json" \
      "$GAPI/commits/$SHA/statuses" \
      2>&1 \
      > "$LOG/read.log"
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
      STATUS="Deploy / $(basename "$VENV") @ $(echo "$HOSTNAME" | grep -o '^[a-z]*')"
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
      echo "{\"state\":\"$STATE\",\"description\":\"Deployed on $TODAY in $TMESS\",\"context\":\"$STATUS\"}" \
      | tee "$LOG/data.json" \
      | curl -s \
         -X POST \
         -H "$AUTH" \
         -D "$LOG/write.txt" \
         -o "$LOG/write.json" \
         "$GAPI/commits/$SHA/statuses" \
         -d @- \
         2>&1 \
         > "$LOG/write.log"

      } &

      wait # DEBUG

   done

   wait
   popd &> /dev/null || exit 12

done \
2>&1 \
| tee "$DISPY_LOG/all.log"

