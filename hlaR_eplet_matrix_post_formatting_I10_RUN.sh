#!/bin/bash
#SBATCH --job-name=hlaR_eplet_matrix_post_formating_I10          # Job Name
#SBATCH --output=/lustre/project/lgragert/gwager/logs/hlaR_eplet_matrix_post_formating_I10_%a.out       ### File in which to store job output
#SBATCH --error=/lustre/project/lgragert/gwager/logs/hlaR_eplet_matrix_post_formating_I10_%a.err 
#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)
#SBATCH --time=24:00:00         # WallTime
#SBATCH --mem=1600
#SBATCH --nodes=1               # Number of Nodes
#SBATCH --ntasks-per-node=1     # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=1       # Number of threads per task (OMP threads)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gwager@tulane.edu
#SBATCH --array=1-14

cd /lustre/project/lgragert/gwager/kidney-outcomes-sfvt/
export CONDA_ENVS_PATH="/lustre/project/lgragert/gwager/conda-envs/"
source activate gwR
unset PYTHONPATH

imp=$(($(($(($SLURM_ARRAY_TASK_ID)) - 1)) / 14))
imp=$(($imp+1))

echo $imp
if [ $imp -eq 1 ]; then
  IMP_ID=1
elif [ $imp -eq 2 ]; then
  IMP_ID=2
elif [ $imp -eq 3 ]; then
  IMP_ID=3
elif [ $imp -eq 4 ]; then
  IMP_ID=4
elif [ $imp -eq 5 ]; then
  IMP_ID=5
elif [ $imp -eq 6 ]; then
  IMP_ID=6
elif [ $imp -eq 7 ]; then
  IMP_ID=7
elif [ $imp -eq 8 ]; then
  IMP_ID=8
elif [ $imp -eq 9 ]; then
  IMP_ID=9
elif [ $imp -eq 10 ]; then
  IMP_ID=10
else
  echo "Fail"
fi 

echo $SLURM_ARRAY_TASK_ID, $IMP_ID
Rscript hlaR_eplet_matrix_post_formating_I10.R $IMP_ID 
