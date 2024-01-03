#!/bin/bash
#SBATCH --job-name=AAMM_RM   ### Job Name
#SBATCH --output=./logs/AAMM_RM_%a.out       ### File in which to store job output
#SBATCH --error=./logs/AAMM_RM_%a.err        ### File in which to store job error messages
#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)
#SBATCH --partition=centos7   ### Required for Python 3.11
#SBATCH --mem=128000
#SBATCH --array=1-20
#SBATCH --time=0-24:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node

replicate=$(($(($SLURM_ARRAY_TASK_ID)) % 10))
if [ $replicate -eq 0 ]; then
  replicate=10
fi

runmatch=$(($(($(($SLURM_ARRAY_TASK_ID)) - 1)) / 10))
matrix=1
if [ $runmatch -eq 1 ]; then
  matrix=0
fi
sfvt=0

echo $SLURM_ARRAY_TASK_ID $replicate $runmatch $matrix $sfvt

SECONDS=0

module unload anaconda3/2020.07
module load anaconda3/2023.07
export CONDA_ENVS_PATH="/lustre/project/lgragert/local/conda-envs/"
source activate lgR
unset PYTHONPATH

echo "Environment set - Executing Python script"

python3 aa_mm_biopython_runmatch_genie_9loc.py $replicate $runmatch $matrix $sfvt

ELAPSED="Elapsed US: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED