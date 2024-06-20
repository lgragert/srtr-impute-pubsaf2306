#!/bin/bash
#BSUB -J AAMM_RM[1-20]                  ### Job Name with array specficied
#BSUB -oo ./logs/AAMM_RM_%I.out         ### File in which to store job output
#BSUB -eo ./logs/AAMM_RM_%I.err         ### File in which to store job error messages
#BSUB -q i2c2_normal                    ### Quality of Service
#BSUB -M 128000                         ### Memory limit per-process (KB)
#BSUB -n 1                              ### sets the core

replicate=$(($((${LSB_JOBINDEX})) % 10))
if [ $replicate -eq 0 ]; then
  replicate=10
fi

runmatch=$(($(($((${LSB_JOBINDEX})) - 1)) / 10))
matrix=1
if [ $runmatch -eq 1 ]; then
  matrix=0
fi
sfvt=0

echo ${LSB_JOBINDEX} $replicate $runmatch $matrix $sfvt

SECONDS=0

module add python/3.7
cd /project/kamoun_shared/code_shared/srtr-impute-pubsaf2306
export CONDA_ENVS_PATH="/project/kamoun_shared/apaynter/local/conda-envs/"
source activate apR
unset PYTHONPATH

echo "Environment set - Executing Python script"

python3.7 aa_mm_biopython_runmatch_genie_9loc.py $replicate $runmatch $matrix $sfvt

ELAPSED="Elapsed US: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED