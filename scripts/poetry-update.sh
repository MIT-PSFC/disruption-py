#!/bin/bash

set -e

PYPROJECT="$(dirname "$0")/../pyproject.toml"
TMPF=$(mktemp)

R="\e[31m"
Y="\e[33m"
G="\e[32m"
B="\e[34m"
U="\e[4m"
X="\e[0m"

{
   echo main
   grep 'tool.poetry.group.[a-z]*.dependencies' "$PYPROJECT" \
   | cut -d . -f 4
} \
| while read -r GROUP
do

   echo -e "\n* $U${GROUP^^}$X *"

   poetry show --latest --top-level --only "$GROUP" \
   | while read -r DEP CURR LAST _
   do

     HAVE=${CURR%.*}
     WANT=${LAST%.*}

     if [[ "$CURR" = "$LAST" ]]
     then

        echo -e "   - $CURR\t=\t$LAST\t$G$DEP$X"

     elif [[ "${HAVE%.*}" != "${WANT%.*}" ]]
     then

        echo -e "   - $R$CURR\t#\t$LAST\t$DEP$X"

     elif [[ "$HAVE" = "$WANT" ]]
     then

        echo -e "   - $CURR\t~\t$LAST\t$B$DEP$X"

     else

        echo -e "   - $Y$CURR\t>\t$LAST\t$DEP$X"
        echo -e "$GROUP\n$DEP@^$WANT.0" >> "$TMPF"

     fi

   done

done

[[ "$1" = "-n" ]] && exit 0
[[ "$1" = "--dry-run" ]] && exit 0

xargs -a "$TMPF" -n 2 poetry add --group
poetry update

git status
git diff --exit-code "$PYPROJECT"
