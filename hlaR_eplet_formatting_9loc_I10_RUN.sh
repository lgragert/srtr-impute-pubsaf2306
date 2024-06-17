#!/bin/bash
#SBATCH --job-name=hlaR_eplet_formatting_9loc_I10                                         ### Job Name
#SBATCH --output=/lustre/project/lgragert/jk/logs/hlaR_eplet_formatting_9loc_%a.out       ### File in which to store job output (based on array task)
#SBATCH --error=/lustre/project/lgragert/jk/logs/hlaR_eplet_formatting_9loc_%a.err        ### File in which to store job error messages (based on array task)
#SBATCH --partition=defq                                                                  ### Partition (default is 'defq')
#SBATCH --qos=normal                                                                      ### Quality of Service (like a queue in PBS)
#SBATCH --mem=1600                                                                        ### Memory to request
#SBATCH --time=0-24:00:00                                                                 ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                                                                         ### Node count required for the job
#SBATCH --ntasks-per-node=1                                                               ### Number of tasks to be launched per Node (MPI processes)
#SBATCH --cpus-per-task=1                                                                 ### Number of threads per task (OMP threads)
#SBATCH --array=1-10                                                                      ### How many jobs to run
#SBATCH --mail-type=ALL                                                                   ### Mails the user the BEGIN, END, FAIL, REQUEUE, or ALL
#SBATCH --mail-user=jko2@tulane.edu                                                       ### tulaneID@tulane.edu

module load R/3.2.4
cd /lustre/project/lgragert/jk/srtr-impute-pubsaf2306/

echo "$SLURM_ARRAY_TASK_ID"                                                               ### SLURM_ARRAY_TASK_ID refers to array number
ARRAY_ID=$(($SLURM_ARRAY_TASK_ID))

Rscript hlaR_eplet_formatting_9loc_I10.R $ARRAY_ID
