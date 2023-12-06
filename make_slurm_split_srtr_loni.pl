#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:  make_slurm_split_srtr_loni.pl
# DESCRIPTION:  Make slurm scripts for running EM pipeline on supercomputer
#
#
# DATE WRITTEN: November 25, 2023
# WRITTEN BY:   Loren Gragert
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ----------    --------------  -------------------------------------
#
##############################################################################
use strict; # always
use warnings; # always
use POSIX;


my $slurm_file = "./slurm/run_slurm_srtr_split_loni.sh";
open (SLURMIMPUTE, ">$slurm_file") || die "Missing $slurm_file\n";


print SLURMIMPUTE "#!/bin/bash\n";
print SLURMIMPUTE "#SBATCH --job-name=SRTR_SPLIT   ### Job Name\n";
print SLURMIMPUTE "#SBATCH --output=./logs/srtr.split.out       ### File in which to store job output\n";
print SLURMIMPUTE "#SBATCH --error=./logs/srtr.split.err        ### File in which to store job error messages\n";
print SLURMIMPUTE "#SBATCH --account=loni_hla_em03   ### Allocation\n";
print SLURMIMPUTE "#SBATCH --partition=single    # Partition default=single\n";
# print SLURMIMPUTE "#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)\n";
print SLURMIMPUTE "#SBATCH --time=0-24:00:00     ### Wall clock time limit in Days-HH:MM:SS\n";
# print SLURMIMPUTE "#SBATCH --mem=128000          ### Memory in MB (default memory is 64000 - 64GB)\n";
print SLURMIMPUTE "#SBATCH --nodes=1             ### Node count required for the job\n";
print SLURMIMPUTE "#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node\n\n\n";

print SLURMIMPUTE "SECONDS=0 # \n\n";

print SLURMIMPUTE "python3 split_pull_pops_srtr.py\n";


print SLURMIMPUTE "ELAPSED=\"Elapsed : \$((\$SECONDS / 3600))hrs \$(((\$SECONDS / 60) % 60))min \$((\$SECONDS % 60))sec\"\n";
print SLURMIMPUTE "echo \$ELAPSED \n\n\n";

close SLURMIMPUTE;
