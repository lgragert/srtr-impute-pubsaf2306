#!/bin/bash
#SBATCH --job-name=AAMM_freq   ### Job Name
#SBATCH --output=./AAMM_freq.out       ### File in which to store job output
#SBATCH --error=./AAMM_freq.err        ### File in which to store job error messages
#SBATCH --qos=long          ### Quality of Service (like a queue in PBS)
#SBATCH --partition=centos7   ### Required for Python 3.11
#SBATCH --mem=256000          ### Memory required because of the DataFrame
#SBATCH --time=0-168:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node
#SBATCH --mail-type=END
#SBATCH --mail-user=apaynter@tulane.edu

module load anaconda3/2023.07
python3 aamm_average_stdev_reps.py
