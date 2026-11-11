#!/bin/bash
set -euo pipefail

module purge
module load openfoam/2312-openmpi5.0.8-gcc14.2.0

script_dir=$(cd "$(dirname "$0")" && pwd)
suite_dir=${1:-/users/sandyz/fastscratch/ichiBoundaryFoam-wave-matrix-261106}
tutorial_dir="$FOAM_TUTORIALS/multiphase/interFoam/laminar/waves"
geometry="$script_dir/../../examples/constant/ichiBoundaryGeometry"

mkdir -p "$suite_dir/cases" "$suite_dir/logs"
: > "$suite_dir/cases.list"

for wave in stokesI stokesII stokesV
do
    case "$wave" in
        stokesI|stokesII)
            end_time=33
            write_interval=3
            ;;
        stokesV)
            end_time=24
            write_interval=2
            ;;
    esac

    for model in laminar kEpsilon kOmegaSST
    do
        for solver in interFoam ichiBoundaryFoam
        do
            name="${wave}-${model}-${solver}"
            case_dir="$suite_dir/cases/$name"
            rm -rf "$case_dir"
            cp -a "$tutorial_dir/$wave" "$case_dir"
            cp -a "$case_dir/0.orig" "$case_dir/0"
            cp "$geometry" "$case_dir/constant/ichiBoundaryGeometry"

            foamDictionary "$case_dir/system/controlDict" \
                -entry application -set "$solver"
            foamDictionary "$case_dir/system/controlDict" \
                -entry endTime -set "$end_time"
            foamDictionary "$case_dir/system/controlDict" \
                -entry writeControl -set adjustableRunTime
            foamDictionary "$case_dir/system/controlDict" \
                -entry writeInterval -set "$write_interval"
            foamDictionary "$case_dir/system/controlDict" \
                -entry purgeWrite -set 0
            foamDictionary "$case_dir/system/controlDict" \
                -entry writePrecision -set 12

            foamDictionary "$case_dir/system/fvSchemes" \
                -entry 'divSchemes.div(rhoPhi,U)' \
                -set 'Gauss limitedLinearV 1'

            if [ "$model" = laminar ]
            then
                foamDictionary \
                    "$case_dir/constant/turbulenceProperties" \
                    -entry simulationType -set laminar
            else
                cp "$script_dir/templates/turbulenceProperties.$model" \
                    "$case_dir/constant/turbulenceProperties"
                cp "$script_dir/templates/k" "$case_dir/0/k"
                cp "$script_dir/templates/nut" "$case_dir/0/nut"

                foamDictionary "$case_dir/system/fvSchemes" \
                    -entry 'divSchemes.div(phi,k)' \
                    -add 'Gauss upwind'
                if [ "$model" = kEpsilon ]
                then
                    cp "$script_dir/templates/epsilon" \
                        "$case_dir/0/epsilon"
                    foamDictionary "$case_dir/system/fvSchemes" \
                        -entry 'divSchemes.div(phi,epsilon)' \
                        -add 'Gauss upwind'
                else
                    cp "$script_dir/templates/omega" \
                        "$case_dir/0/omega"
                    foamDictionary "$case_dir/system/fvSchemes" \
                        -entry 'divSchemes.div(phi,omega)' \
                        -add 'Gauss upwind'
                    foamDictionary "$case_dir/system/fvSchemes" \
                        -entry wallDist \
                        -add '{ method meshWave; }'
                fi

                foamDictionary "$case_dir/system/fvSolution" \
                    -entry 'solvers."(k|epsilon|omega)"' \
                    -add '{ solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.1; }'
                foamDictionary "$case_dir/system/fvSolution" \
                    -entry 'solvers."(k|epsilon|omega)Final"' \
                    -add '{ solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0; }'
            fi

            printf '%s\n' "$name" >> "$suite_dir/cases.list"
        done
    done
done

cp "$script_dir/run-array.sbatch" "$suite_dir/run-array.sbatch"
cp "$script_dir/compare.py" "$suite_dir/compare.py"
printf 'Prepared %s cases in %s\n' \
    "$(wc -l < "$suite_dir/cases.list")" "$suite_dir"
