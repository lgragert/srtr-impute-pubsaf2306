#!/bin/bash
#SBATCH --job-name=hlaR_eplet_formating_9loc_I10          # Job Name
#SBATCH --output=/lustre/project/lgragert/gwager/logs/hlaR_eplet_formating_9loc_%a.out       ### File in which to store job output
#SBATCH --error=/lustre/project/lgragert/gwager/logs/hlaR_eplet_formating_9loc_%a.err 
#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)
#SBATCH --time=24:00:00         # WallTime
#SBATCH --mem=1600
#SBATCH --nodes=1               # Number of Nodes
#SBATCH --ntasks-per-node=1     # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=1       # Number of threads per task (OMP threads)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gwager@tulane.edu
#SBATCH --array=1-10
cd /lustre/project/lgragert/gwager/kidney-outcomes-sfvt/
export CONDA_ENVS_PATH="/lustre/project/lgragert/gwager/conda-envs/"
source activate gwR
unset PYTHONPATH

echo "$SLURM_ARRAY_TASK_ID"
ARRAY_ID=$(($SLURM_ARRAY_TASK_ID))

Rscript hlaR_eplet_formating_9loc_I10.R $ARRAY_ID