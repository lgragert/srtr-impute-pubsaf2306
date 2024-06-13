#!/bin/bash
#BSUB -J AAMM_freq   ### Job Name
#BSUB -o /project/kamoun_shared/apaynter/logs/AAMM_freq.out       ### File in which to store job output
#BSUB -e /project/kamoun_shared/apaynter/logs/AAMM_freq.err        ### File in which to store job error messages
#BSUB -q i2c2_normal          ### Quality of Service
#BSUB -M 256000          ### Memory required
#BSUB -n 1             ### Node count required for the job
#BSUB -u apaynter@tulane.edu
#BSUB -N		### Emails when job ends

module add python/3.7
cd /project/kamoun_shared/apaynter/
python3.7 aamm_average_stdev_reps.py
