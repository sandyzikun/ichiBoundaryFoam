#!/bin/bash -l
#SBATCH --job-name=ichiWaveMatrix
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --array=0-17%6
#SBATCH --output=logs/%A-%a.out
#SBATCH --error=logs/%A-%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?}"
module purge
module load openfoam/2312-openmpi5.0.8-gcc14.2.0
export PATH="/users/sandyz/slurm_of/waves2-3d/ichiBoundaryFoam/platforms/linux64GccDPInt32-spack/bin:$PATH"
export FOAM_SIGFPE=true
export FOAM_SETNAN=true

name=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" cases.list)
case_dir="$SLURM_SUBMIT_DIR/cases/$name"
application=$(foamDictionary "$case_dir/system/controlDict" \
    -entry application -value)

cd "$case_dir"
blockMesh > log.blockMesh 2>&1
setFields > log.setFields 2>&1
"$application" > "log.$application" 2>&1

foamListTimes -latestTime > latestTime.txt
touch RUN_COMPLETE
