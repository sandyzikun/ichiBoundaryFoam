#!/usr/bin/env bash
set -euo pipefail
##
suite_dir="${1:-/users/sandyz/fastscratch/ichiBoundaryFoam-validation-261104}"
printf '%-24s %12s %12s %12s %12s\n' \
    case final_time slip div phase_drift load_x
##
while IFS= read -r name; do
    log="$suite_dir/cases/$name/log.run"
    if [[ ! -f "$log" ]]; then
        printf '%-24s %12s\n' "$name" MISSING
        continue
    fi
    final_time=$(awk '/^Time = /{v=$3} END{print v+0}' "$log")
    slip=$(awk \
        '/IBM projection residuals:/{v=$7} END{print v+0}' \
        "$log")
    div=$(awk \
        '/IBM projection residuals:/{v=$11} END{print v+0}' \
        "$log")
    phase_drift=$(awk '
        /IBM accessible phase volume:/{
            v=$NF+0; if (!seen++) first=v; last=v
        }
        END{
            if (seen && first != 0) print (last-first)/first;
            else print 0
        }' "$log")
    load_x=$(awk '
        /IBM preliminary load:/{
            line=$0; sub(/^.*force = \(/, "", line);
            split(line, a, " "); v=a[1]
        }
        END{print v+0}' "$log")
    printf '%-24s %12g %12g %12g %12g %12g\n' \
        "$name" "$final_time" "$slip" "$div" \
        "$phase_drift" "$load_x"
done < "$suite_dir/case-list.txt"

