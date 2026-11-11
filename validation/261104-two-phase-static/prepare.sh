#!/usr/bin/env bash
set -euo pipefail
##
project_dir="/users/sandyz/slurm_of/waves2-3d/ichiBoundaryFoam"
template="$project_dir/validation/261103-cylinder/cases/ichi-ibm-stokesII-mdf"
suite_dir="${1:-/users/sandyz/fastscratch/ichiBoundaryFoam-2phase-static-261104}"
solver="$project_dir/platforms/linux64GccDPInt32-spack/bin/ichiBoundaryFoam"
##
if [[ -e "$suite_dir/PREPARED_AT" ]]; then
    echo "Suite already exists: $suite_dir" >&2
    exit 1
fi
if [[ ! -x "$solver" ]]; then
    echo "Solver not found: $solver" >&2
    exit 1
fi
##
mkdir -p "$suite_dir/cases"
cp "$solver" "$suite_dir/ichiBoundaryFoam"
cp "$project_dir/validation/261104-two-phase-static/sbatch.sh" \
    "$suite_dir/"
cp "$project_dir/validation/261104-two-phase-static/summarise.sh" \
    "$suite_dir/"
cp "$project_dir/validation/261104-two-phase-static/README.md" \
    "$suite_dir/"
##
copy_case()
{
    local name="$1"
    local enabled="$2"
    local case_dir="$suite_dir/cases/$name"
    ##
    mkdir -p "$case_dir"
    cp -R "$template/0.orig" "$case_dir/"
    cp -R "$template/constant" "$case_dir/"
    cp -R "$template/system" "$case_dir/"
    rm -rf "$case_dir/constant/polyMesh"
    ##
    sed -i -E \
        's/^startFrom[[:space:]]+.*/startFrom       startTime;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^endTime[[:space:]]+.*/endTime         9;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^deltaT[[:space:]]+.*/deltaT          0.0025;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^writeInterval[[:space:]]+.*/writeInterval   0.25;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^maxCo[[:space:]]+.*/maxCo           0.4;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^maxAlphaCo[[:space:]]+.*/maxAlphaCo      0.4;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^maxDeltaT[[:space:]]+.*/maxDeltaT       0.01;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        "s/enabled[[:space:]]+(true|false);/enabled         $enabled;/" \
        "$case_dir/constant/ichiBoundaryGeometry"
    sed -i -E \
        's/(nMultiDirectForcingCorrectors).*/\1 4;/' \
        "$case_dir/constant/ichiBoundaryGeometry"
    sed -i -E \
        's/(slipTolerance).*/\1 1e-5;/' \
        "$case_dir/constant/ichiBoundaryGeometry"
}
##
copy_case stokesII-control false
copy_case stokesII-static-ibm true
printf '%s\n' stokesII-control stokesII-static-ibm \
    > "$suite_dir/case-list.txt"
date --iso-8601=seconds > "$suite_dir/PREPARED_AT"
echo "Prepared: $suite_dir"

