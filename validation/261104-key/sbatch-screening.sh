#!/bin/bash -l
#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH -J ichiIbmKey
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-7
#SBATCH -t 0-06:00:00
set -euo pipefail
##
cd "${SLURM_SUBMIT_DIR:?'variable not set'}"
ulimit -s unlimited
export OMP_NUM_THREADS=1
module purge
module load openfoam/2312-openmpi5.0.8-gcc14.2.0
##
mapfile -t cases < case-list.txt
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
case_dir="$SLURM_SUBMIT_DIR/cases/$case_name"
solver="$SLURM_SUBMIT_DIR/ichiBoundaryFoam"
##
cd "$case_dir"
echo "Case: $case_name"
echo "Started: $(date --iso-8601=seconds)"
if [[ -d constant/polyMesh ]]; then
    echo "Using the prepared mesh" > log.blockMesh
else
    blockMesh > log.blockMesh 2>&1
fi
checkMesh > log.checkMesh 2>&1
if [[ "$case_name" == body-fitted-* ]]; then
    "$solver" > log.run 2>&1
else
    "$solver" > log.run 2>&1
fi
echo "Finished: $(date --iso-8601=seconds)"

