#!/bin/bash
#SBATCH --job-name=TRS_AA   ### Job Name
#SBATCH --output=./logs/TRS_AA_%a.out       ### File in which to store job output
#SBATCH --error=./logs/TRS_AA_%a.err        ### File in which to store job error messages
#SBATCH --qos=long          ### Quality of Service (like a queue in PBS)
#SBATCH --partition=centos7   ### Required for Python 3.11
#SBATCH --array=1-6
#SBATCH --mem=128000
#SBATCH --time=0-168:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node

module unload anaconda3/2020.07
module load anaconda3/2023.07
export CONDA_ENVS_PATH="/lustre/project/lgragert/local/conda-envs/"
conda activate lgR

pop_id=$(($(($(($SLURM_ARRAY_TASK_ID)) - 1)) / 6))
pop_id=$(($pop_id+1))

echo $pop_id
if [ $pop_id -eq 1 ]; then
  pop="AFA"
elif [ $pop_id -eq 2 ]; then
  pop="API"
elif [ $pop_id -eq 3 ]; then
  pop="CAU"
elif [ $pop_id -eq 4 ]; then
  pop="HIS"
elif [ $pop_id -eq 5 ]; then
  pop="NAM"
elif [ $pop_id -eq 6 ]; then
  pop="MLT"
else
  echo "Pop out of range"
fi 

python3 typing_res_score_9loc_AA.py $pop
