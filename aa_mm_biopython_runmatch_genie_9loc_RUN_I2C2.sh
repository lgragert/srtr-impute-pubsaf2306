#!/bin/bash
#BSUB -J AAMM_RM                 ### Job Name with array specficied
#BSUB -o ./logs/%J.out         ### File in which to store job output
#BSUB -e ./logs/%J.err         ### File in which to store job error messages
#BSUB -q i2c2_normal                    ### Quality of Service
#BSUB -M 128000                         ### Memory limit per-process
#BSUB -u apaynter@tulane.edu            ### send me mail
#BSUB -N

source activate /project/kamoun_shared/apaynter/conda-envs/apR
cd /project/kamoun_shared/code_shared/srtr-impute-pubsaf2306

python aa_mm_biopython_runmatch_genie_9loc.py 1 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 1 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 2 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 2 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 3 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 3 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 4 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 4 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 5 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 5 0 1 0

python aa_mm_biopython_runmatch_genie_9loc.py 6 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 6 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 7 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 7 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 8 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 8 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 9 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 9 0 1 0
python aa_mm_biopython_runmatch_genie_9loc.py 10 1 0 0
python aa_mm_biopython_runmatch_genie_9loc.py 10 0 1 0

ELAPSED="Elapsed US: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
