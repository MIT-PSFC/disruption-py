#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# disruption-errors.sh
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# this script parses a `disruption-py` log file and collects unexpected errors.
#
# first, it runs `disruption-py` by passing on any command-line argument, but
# can optionally digest a previous log file, instead.
# then, it parses the log file and collects all errors and critical messages,
# excluding those that match the patterns in the sibling .ignore file.
# finally, it prints a list of errors, together with their relative frequency.
#
# usage examples:
#
#    disruption-errors.sh -m get_ip_parameters disruption_warning
#    disruption-errors.sh /path/to/run/output.log
#    ls /path/to/run/*log | xargs -n1 disruption-errors.sh
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

IGNORE="${BASH_SOURCE[0]%.sh}.ignore"

TMPD="${LOCALSCRATCH:-/tmp}/$USER/disruption-py/.$(date +%F)"
mkdir -p "$TMPD" || exit 10
TMPF=$(mktemp -p "$TMPD" "errors-$(date +%s)-XXX.log")
TMPS=${TMPF%.log}.txt

if [[ $# -eq 1 ]] && [[ -f "$1" ]] && [[ "$1" =~ \.log$ ]]
then

   LOG=$1

else

   uv run disruption-py -l info "$@" \
   | grep -e ERROR -e CRITICAL -e 'INFO.*workflow' -e 'INFO.*Logging' \
   | grep -vFf "$IGNORE" \
   | tee "$TMPF"
   echo "disruption-py: rc ${PIPESTATUS[0]}"

   LOG=$(grep -o 'Logging:.*' "$TMPF" | cut -d' ' -f2)

fi

[[ -n "$LOG" ]] || exit 11
[[ -s "$LOG" ]] || exit 12

echo -e "\033[36m"
realpath -m "$LOG" "$TMPF" "$TMPS"

echo -e "\033[31m"
grep -e ERROR -e CRITICAL "$LOG" \
| grep -vFf "$IGNORE" \
| cut -d'|' -f2 \
| sort \
| uniq -c \
| sort -n \
| tee "$TMPS"
echo -e "\033[0m"

[[ ! -s "$TMPS" ]]
