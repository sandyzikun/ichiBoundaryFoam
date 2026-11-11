#!/bin/bash -l
#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH -J ichi2pStatic
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1
#SBATCH -t 0-04:00:00
set -euo pipefail
##
cd "${SLURM_SUBMIT_DIR:?'variable not set'}"
ulimit -s unlimited
export OMP_NUM_THREADS=1
module purge
module load openfoam/2312-openmpi5.0.8-gcc14.2.0
##
mapfile -t cases < case-list.txt
name="${cases[$SLURM_ARRAY_TASK_ID]}"
case_dir="$SLURM_SUBMIT_DIR/cases/$name"
solver="$SLURM_SUBMIT_DIR/ichiBoundaryFoam"
##
cd "$case_dir"
echo "Case: $name"
echo "Started: $(date --iso-8601=seconds)"
cp -R 0.orig 0
blockMesh > log.blockMesh 2>&1
setFields > log.setFields 2>&1
"$solver" > log.run 2>&1
echo "Finished: $(date --iso-8601=seconds)"

