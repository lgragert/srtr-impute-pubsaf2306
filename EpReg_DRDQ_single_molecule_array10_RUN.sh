#!/bin/bash
#SBATCH --job-name=EpReg_DRDQ_single_molecule_array10                    ### Job Name
#SBATCH --array=1-10                                                     ### Number of array tasks
#SBATCH --output=./logs/EpReg_DRDQ_single_molecule_array10_%A_%a.out     ### File in which to store job output (based on array task)
#SBATCH --error=./logs/EpReg_DRDQ_single_molecule_array10_%A_%a.err      ### File in which to store job error messages (based on array task)
#SBATCH --partition=centos7                                              ### Partition required for Python 3.11 (default is 'defq')
#SBATCH --qos=normal                                                     ### Quality of Service (like a queue in PBS)
#SBATCH --time=0-24:00:00                                                ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --mem=256000                                                     ### Memory required because of the DataFrame




#SBATCH --mail-type=ALL                                                  ### Mails the user the BEGIN, END, FAIL, REQUEUE, or ALL
#SBATCH --mail-user=jko2@tulane.edu                                      ### tulaneID@tulane.edu

module load anaconda3/2023.07
cd /lustre/project/lgragert/jk/srtr-impute-pubsaf2306/

echo "Running task with SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"       ### SLURM_ARRAY_TASK_ID refers to array number

python3 EpReg_DRDQ_single_molecule_array10.py
