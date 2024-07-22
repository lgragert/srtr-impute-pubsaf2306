#!/bin/bash
#SBATCH --job-name=hla_epitope_registry_eplet_test                                        ### Job Name
#SBATCH --output=/lustre/project/lgragert/jk/logs/hla_epitope_registry_eplet_test.out     ### File in which to store job output (based on array task)
#SBATCH --error=/lustre/project/lgragert/jk/logs/hla_epitope_registry_eplet_test.err      ### File in which to store job error messages (based on array task)
#SBATCH --partition=defq                                                                  ### Partition (default is 'defq')
#SBATCH --qos=long                                                                        ### Quality of Service (like a queue in PBS)
#SBATCH --mem=16000                                                                       ### Memory to request
#SBATCH --time=0-48:00:00                                                                 ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                                                                         ### Node count required for the job
#SBATCH --ntasks-per-node=1                                                               ### Number of tasks to be launched per Node (MPI processes)
#SBATCH --cpus-per-task=1                                                                 ### Number of threads per task (OMP threads)

#SBATCH --mail-type=ALL                                                                   ### Mails the user the BEGIN, END, FAIL, REQUEUE, or ALL
#SBATCH --mail-user=jko2@tulane.edu                                                       ### tulaneID@tulane.edu

module load anaconda3/2023.07
cd /lustre/project/lgragert/jk/srtr-impute-pubsaf2306/
python3 hla_epitope_registry_eplet_test.py
