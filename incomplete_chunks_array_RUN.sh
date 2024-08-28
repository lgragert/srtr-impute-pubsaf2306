#!/bin/bash
#SBATCH --job-name=incomplete_chunks_array                               ### Job Name
#SBATCH --array=1-10                                                     ### Number of array tasks
#SBATCH --output=./logs/incomplete_chunks_array_%A_%a.out                ### File in which to store job output (based on array task)
#SBATCH --error=./logs/incomplete_chunks_array_%A_%a.err                 ### File in which to store job error messages (based on array task)
#SBATCH --partition=centos7                                              ### Partition required for Python 3.11 (default is 'defq')
#SBATCH --qos=long                                                       ### Quality of Service (like a queue in PBS)
#SBATCH --mem=256000                                                     ### Memory required because of the DataFrame
#SBATCH --time=0-168:00:00                                               ### Wall clock time limit in Days-HH:MM:SS




#SBATCH --mail-type=ALL                                                  ### Mails the user the BEGIN, END, FAIL, REQUEUE, or ALL
#SBATCH --mail-user=jko2@tulane.edu                                      ### tulaneID@tulane.edu

module load anaconda3/2023.07
cd /lustre/project/lgragert/jk/srtr-impute-pubsaf2306/

echo "Running task with SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"       ### SLURM_ARRAY_TASK_ID refers to array number

python3 incomplete_chunks_array.py
