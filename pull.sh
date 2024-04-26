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

# API token
GAPI="https://api.github.com/repos/mit-psfc/disruption-py"
AUTH="$(cat "/home/$USER/.gh_pat")"

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

   # read statuses
   SHA=$(git rev-parse HEAD)
   if [[ -n "$SHA" ]]
   then
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

      {

      # read status
      STATUS="Install / py${VENV##*py} @ ${HOSTNAME%-*}"
      if [[ -s "$LOG/sha.json" ]]
      then
         STATE=$(jq -r ".[]|select(.context==\"$STATUS\").state" "$LOG/sha.json")
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

      # run
      if [[ -s examples/efit.py ]]
      then
         poetry run python examples/efit.py \
         1> "$LOG/efit.out" \
         2> "$LOG/efit.err"
      elif [[ -s examples/quick.py ]]
      then
         poetry run python examples/quick.py \
         1> "$LOG/quick.out" \
         2> "$LOG/quick.err"
      else
         poetry run python -c "import disruption_py as dpy; print(dpy.__file__)" \
         1> "$LOG/import.out" \
         2> "$LOG/import.err"
      fi

      # rc
      RC=${PIPESTATUS[0]}
      [[ $RC -eq 0 ]] && STATE=success || STATE=failure

      # log
      echo -e "$(date) :: $DISPY_BRANCH @ py$DISPY_PYVERS = rc $RC"

      # status
      echo "{\"state\":\"$STATE\",\"description\":\"Updated in ${SECONDS} s\",\"context\":\"$STATUS\"}" \
      | tee "$LOG/data.json" \
      | curl -s \
         -X POST \
         -H "$AUTH" \
         -D "$LOG/write.txt" \
         -o "$LOG/write.json" \
         "$GAPI/commits/$SHA/statuses" \
         -d @- \
         2>&1 \
         >> "$LOG/git.log"

      # deactivate
      deactivate

      } &

   done

   wait
   popd &> /dev/null || exit 12

done \
2>&1 \
| tee "$DISPY_LOG/all.log"

