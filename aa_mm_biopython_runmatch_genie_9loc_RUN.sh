#!/bin/bash
#SBATCH --job-name=AAMM_RM   ### Job Name
#SBATCH --output=./logs/aa_mm_biopython_runmatch_9loc.out       ### File in which to store job output
#SBATCH --error=./logs/aa_mm_biopython_runmatch_9loc.err        ### File in which to store job error messages
#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)
#SBATCH --partition=centos7   ### Required for Python 3.11
#SBATCH --mem=256000
#SBATCH --time=0-168:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node

module unload anaconda3/2020.07
module unload anaconda3/2023.07
export CONDA_ENVS_PATH="/lustre/project/lgragert/local/conda-envs/"
source activate lgR
unset PYTHONPATH
python3 aa_mm_biopython_runmatch_genie_9loc.py
