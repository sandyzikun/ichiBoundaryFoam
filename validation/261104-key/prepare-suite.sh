#!/usr/bin/env bash
set -euo pipefail
##
project_dir="/users/sandyz/slurm_of/waves2-3d/ichiBoundaryFoam"
source_dir="$project_dir/validation/261103-cylinder"
suite_dir="${1:-/users/sandyz/fastscratch/ichiBoundaryFoam-validation-261104}"
solver="$project_dir/platforms/linux64GccDPInt32-spack/bin/ichiBoundaryFoam"
##
if [[ ! -x "$solver" ]]; then
    echo "Solver not found: $solver" >&2
    exit 1
fi
##
if [[ -e "$suite_dir/PREPARED_AT" ]]; then
    echo "Suite already exists: $suite_dir" >&2
    echo "Move or archive it before preparing a new suite." >&2
    exit 1
fi
##
mkdir -p "$suite_dir/cases"
cp "$project_dir/validation/261104-key/sbatch-screening.sh" \
    "$suite_dir/"
cp "$project_dir/validation/261104-key/summarise.sh" \
    "$suite_dir/"
cp "$project_dir/validation/261104-key/README.md" \
    "$suite_dir/README.md"
cp "$solver" "$suite_dir/ichiBoundaryFoam"
##
copy_ibm_case()
{
    local name="$1"
    local nx="$2"
    local ny="$3"
    local dt="$4"
    local end_time="$5"
    local width="$6"
    local motion="$7"
    local case_dir="$suite_dir/cases/$name"
    ##
    mkdir -p "$case_dir"
    cp -R "$source_dir/immersed/0" "$case_dir/"
    cp -R "$source_dir/immersed/constant" "$case_dir/"
    cp -R "$source_dir/immersed/system" "$case_dir/"
    cp "$source_dir/common/system/fvSchemes" "$case_dir/system/"
    cp "$source_dir/common/system/fvSolution" "$case_dir/system/"
    cp "$source_dir/common/system/probes" "$case_dir/system/"
    ##
    sed -i -E \
        "s/\(173 120 1\)/($nx $ny 1)/" \
        "$case_dir/system/blockMeshDict"
    sed -i -E \
        "s/^endTime[[:space:]]+.*/endTime         $end_time;/" \
        "$case_dir/system/controlDict"
    sed -i -E \
        "s/^deltaT[[:space:]]+.*/deltaT          $dt;/" \
        "$case_dir/system/controlDict"
    sed -i -E \
        "s/^adjustTimeStep[[:space:]]+.*/adjustTimeStep  no;/" \
        "$case_dir/system/controlDict"
    sed -i -E \
        "s/^maxDeltaT[[:space:]]+.*/maxDeltaT       $dt;/" \
        "$case_dir/system/controlDict"
    sed -i -E \
        "s/transitionWidth[[:space:]]+[0-9.eE+-]+;/transitionWidth     $width;/" \
        "$case_dir/constant/ichiBoundaryGeometry"
    sed -i \
        '/transitionWidth/a\    multiDirectForcing  true;' \
        "$case_dir/constant/ichiBoundaryGeometry"
    sed -i \
        '/multiDirectForcing/a\    nMultiDirectForcingCorrectors 4;' \
        "$case_dir/constant/ichiBoundaryGeometry"
    sed -i \
        '/nMultiDirectForcingCorrectors/a\    slipTolerance 1e-5;' \
        "$case_dir/constant/ichiBoundaryGeometry"
    ##
    if [[ "$motion" == "moving" ]]; then
        sed -i \
            's/type                static;/type                constant6DoF;/' \
            "$case_dir/constant/ichiBoundaryGeometry"
        sed -i \
            '/referenceTime/a\    translationVelocity (0 -0.002 0);' \
            "$case_dir/constant/ichiBoundaryGeometry"
        sed -i \
            '/translationVelocity/a\    rotationRate (0 0 0);' \
            "$case_dir/constant/ichiBoundaryGeometry"
    fi
}
##
copy_body_fitted_case()
{
    local case_dir="$suite_dir/cases/body-fitted-medium"
    ##
    mkdir -p "$case_dir"
    cp -R "$source_dir/body-fitted/0" "$case_dir/"
    cp -R "$source_dir/body-fitted/constant" "$case_dir/"
    cp -R "$source_dir/body-fitted/system" "$case_dir/"
    cp "$source_dir/common/system/fvSchemes" "$case_dir/system/"
    cp "$source_dir/common/system/fvSolution" "$case_dir/system/"
    cp "$source_dir/common/system/probes" "$case_dir/system/"
    cp "$source_dir/common/system/forcesBody" "$case_dir/system/"
    cp -R \
        "$source_dir/cases/ichi-body-fitted/constant/polyMesh" \
        "$case_dir/constant/"
    ##
    sed -i -E \
        's/^endTime[[:space:]]+.*/endTime         0.5;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^deltaT[[:space:]]+.*/deltaT          0.00125;/' \
        "$case_dir/system/controlDict"
    sed -i -E \
        's/^adjustTimeStep[[:space:]]+.*/adjustTimeStep  no;/' \
        "$case_dir/system/controlDict"
    sed -i '$a\functions' "$case_dir/system/controlDict"
    sed -i '$a\{' "$case_dir/system/controlDict"
    sed -i '$a\    #include "forcesBody"' \
        "$case_dir/system/controlDict"
    sed -i '$a\    #include "probes"' \
        "$case_dir/system/controlDict"
    sed -i '$a\}' "$case_dir/system/controlDict"
}
##
copy_ibm_case space-coarse 86 60 0.00125 0.5 0.02 static
copy_ibm_case space-medium 173 120 0.00125 0.5 0.01 static
copy_ibm_case space-fine 346 240 0.00125 0.5 0.005 static
copy_ibm_case time-dt005 173 120 0.005 0.5 0.01 static
copy_ibm_case time-dt0025 173 120 0.0025 0.5 0.01 static
copy_ibm_case time-dt00125 173 120 0.00125 0.5 0.01 static
copy_ibm_case moving-medium 173 120 0.0025 5 0.01 moving
copy_body_fitted_case
##
date --iso-8601=seconds > "$suite_dir/PREPARED_AT"
printf '%s\n' \
    space-coarse \
    space-medium \
    space-fine \
    time-dt005 \
    time-dt0025 \
    time-dt00125 \
    body-fitted-medium \
    moving-medium > "$suite_dir/case-list.txt"
echo "Prepared: $suite_dir"

