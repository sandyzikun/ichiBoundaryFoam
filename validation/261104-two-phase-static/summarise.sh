#!/usr/bin/env bash
set -euo pipefail
##
suite_dir="${1:-/users/sandyz/fastscratch/ichiBoundaryFoam-2phase-static-261104}"
printf '%-24s %8s %12s %12s %12s %12s %12s\n' \
    case time min_alpha max_alpha phase_frac slip divergence
##
while IFS= read -r name; do
    log="$suite_dir/cases/$name/log.run"
    if [[ ! -f "$log" ]]; then
        printf '%-24s %8s\n' "$name" MISSING
        continue
    fi
    awk -v name="$name" '
        /^Time = / {time=$3}
        /Phase-1 volume fraction =/ {
            phase=$5; minAlpha=$8; maxAlpha=$11
        }
        /IBM projection residuals:/ {
            slip=$7; divergence=$11
        }
        END {
            printf "%-24s %8g %12g %12g %12g %12g %12g\n", \
                name, time, minAlpha, maxAlpha, phase, slip, divergence
        }
    ' "$log"
done < "$suite_dir/case-list.txt"

