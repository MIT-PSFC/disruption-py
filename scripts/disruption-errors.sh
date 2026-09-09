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
# finally, it prints a useful summary, consisting of:
# - a (red) list of unique errors, together with their relative frequency;
# - a (yellow) list of error examples, one per unique error;
# - a (green) list of shot numbers, one per unique error.
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
TMPS=${TMPF%.log}.stats
TMPE=${TMPF%.log}.errs
TMPL=${TMPF%.log}.shots

if [[ $# -eq 1 ]] && [[ -f "$1" ]] && [[ "$1" =~ \.log$ ]]
then

   LOG=$1

else

   uv run disruption-py -l info "$@" \
   | grep -e ERROR -e CRITICAL -e 'INFO.*workflow' -e 'INFO.*Logging' \
   | grep -vFf "$IGNORE" \
   | tee "$TMPF"

   LOG=$(grep -o 'Logging:.*\.log' "$TMPF" | cut -d' ' -f2)

fi

[[ -n "$LOG" ]] || exit 11
[[ -s "$LOG" ]] || exit 12

echo -e "\033[36m"
realpath -m "$LOG" "$TMPF" "$TMPE" "$TMPS" "$TMPL"

echo -e "\033[31m"
grep -e ERROR -e CRITICAL "$LOG" \
| grep -vFf "$IGNORE" \
| cut -d'|' -f2 \
| sort \
| tee "$TMPE" \
| uniq -c \
| sort -n \
| tee "$TMPS"

echo -e "\033[33m"
sort -u "$TMPE" -o "$TMPE"
while IFS= read -r ERR; do
  grep -Fm1 "$ERR" "$LOG" \
  | sed 's/^[^\[]*//'
done < "$TMPE" \
| tee "$TMPL.2"

echo -e "\033[32m"
cut -d'|' -f1 "$TMPL.2" \
| cut -d'#' -f2 \
| sort -un \
| tee "$TMPL" \
| xargs

echo -e "\033[0m"
[[ ! -s "$TMPS" ]]
