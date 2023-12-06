#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:  make_slurm_EM.pl
# DESCRIPTION:  Make slurm scripts for running EM pipeline on supercomputer
#
#
# DATE WRITTEN: April 19, 2020
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

my $pop_mode = shift @ARGV;

# different sets of populations depending on how much time you have - default is US-only NEMO
my $pops;
$pops = "CAU,HIS,AFA,ASN,NAM,MLT";

my $dataset = "srtr";

my @pops = split /,/,$pops;


my $sbatch_file = "run_sbatch_impute_srtr_loni.sh";
open (SBATCH, ">$sbatch_file") || die "Missing $sbatch_file\n";


foreach my $pop (@pops) {

  my $pop_greedy = $pop;

  my $slurm_impute_file = "./slurm/run_slurm_srtr_" . $pop . "_impute_loni.sh";
  open (SLURMIMPUTE, ">$slurm_impute_file") || die "Missing $slurm_impute_file\n";

  print SBATCH "jid_impute_$pop=\$(sbatch $slurm_impute_file)\n";
  print SBATCH "echo \$jid_impute_$pop impute $pop\n";
  print SBATCH "echo \${jid_impute_$pop\#\#\* }\n";

  my $job_impute = $pop . "_impute";


  print SLURMIMPUTE "#!/bin/bash\n";
  print SLURMIMPUTE "#SBATCH --job-name=$job_impute   ### Job Name\n";
  print SLURMIMPUTE "#SBATCH --output=./logs/srtr_$pop.impute.out       ### File in which to store job output\n";
  print SLURMIMPUTE "#SBATCH --error=./logs/srtr_$pop.impute.err        ### File in which to store job error messages\n";
  print SLURMIMPUTE "#SBATCH --account=loni_hla_em03   ### Allocation\n";
  print SLURMIMPUTE "#SBATCH --partition=single    # Partition default=single\n";
  # print SLURMIMPUTE "#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)\n";
  if ($pop eq "CAU") {
    print SLURMIMPUTE "#SBATCH --time=0-168:00:00     ### Wall clock time limit in Days-HH:MM:SS\n";
  }
  else {
    print SLURMIMPUTE "#SBATCH --time=0-24:00:00     ### Wall clock time limit in Days-HH:MM:SS\n"; 
  }
  if (($pop eq "AFA") || ($pop eq "CAU")) {
    print SLURMIMPUTE "#SBATCH --mem=128000          ### Memory in MB (default memory is 64000 - 64GB)\n";
  }
  print SLURMIMPUTE "#SBATCH --nodes=1             ### Node count required for the job\n";
  print SLURMIMPUTE "#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node\n\n\n";

  print SLURMIMPUTE "SECONDS=0 # $pop\n\n";

  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop C B greedy 0 CB 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./CB/pull.$dataset.$pop.txt ./CB/glid.$dataset.$pop.txt ./CB/freqs.$pop.csv ./CB/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop A C~B greedy 0 ACB 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./ACB/pull.$dataset.$pop.txt ./ACB/glid.$dataset.$pop.txt ./ACB/freqs.$pop.csv ./ACB/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop DRBX DRB1 greedy 0 DRBXDRB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./DRBXDRB1/pull.$dataset.$pop.txt ./DRBXDRB1/glid.$dataset.$pop.txt ./DRBXDRB1/freqs.$pop.csv ./DRBXDRB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop DRBX~DRB1 DQB1 greedy 0 DRBXDRB1DQB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./DRBXDRB1DQB1/pull.$dataset.$pop.txt ./DRBXDRB1DQB1/glid.$dataset.$pop.txt ./DRBXDRB1DQB1/freqs.$pop.csv ./DRBXDRB1DQB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop DRBX~DRB1~DQB1 DQA1 greedy 0 DRBXDRB1DQA1DQB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./DRBXDRB1DQA1DQB1/pull.$dataset.$pop.txt ./DRBXDRB1DQA1DQB1/glid.$dataset.$pop.txt ./DRBXDRB1DQA1DQB1/freqs.$pop.csv ./DRBXDRB1DQA1DQB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop DRBX~DRB1~DQB1 DPB1 greedy 0 DRBXDRB1DQB1DPB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./DRBXDRB1DQB1DPB1/pull.$dataset.$pop.txt ./DRBXDRB1DQB1DPB1/glid.$dataset.$pop.txt ./DRBXDRB1DQB1DPB1/freqs.$pop.csv ./DRBXDRB1DQB1DPB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop DRBX~DRB1~DQA1~DQB1 DPB1 greedy 0 DRBXDRB1DQA1DQB1DPB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./DRBXDRB1DQA1DQB1DPB1/pull.$dataset.$pop.txt ./DRBXDRB1DQA1DQB1DPB1/glid.$dataset.$pop.txt ./DRBXDRB1DQA1DQB1DPB1/freqs.$pop.csv ./DRBXDRB1DQA1DQB1DPB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop DRBX~DRB1~DQA1~DQB1~DPB1 DPA1 greedy 0 DRBXDRB1DQA1DQB1DPA1DPB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./DRBXDRB1DQA1DQB1DPA1DPB1/pull.$dataset.$pop.txt ./DRBXDRB1DQA1DQB1DPA1DPB1/glid.$dataset.$pop.txt ./DRBXDRB1DQA1DQB1DPA1DPB1/freqs.$pop.csv ./DRBXDRB1DQA1DQB1DPA1DPB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop A~C~B DRBX~DRB1~DQA1~DQB1~DPA1~DPB1 greedy 0 ACBDRBXDRB1DQA1DQB1DPA1DPB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./ACBDRBXDRB1DQA1DQB1DPA1DPB1/pull.$dataset.$pop.txt ./ACBDRBXDRB1DQA1DQB1DPA1DPB1/glid.$dataset.$pop.txt ./ACBDRBXDRB1DQA1DQB1DPA1DPB1/freqs.$pop.csv ./ACBDRBXDRB1DQA1DQB1DPA1DPB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop A~C~B DRBX~DRB1~DQB1 greedy 0 ACBDRBXDRB1DQB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./ACBDRBXDRB1DQB1/pull.$dataset.$pop.txt ./ACBDRBXDRB1DQB1/glid.$dataset.$pop.txt ./ACBDRBXDRB1DQB1/freqs.$pop.csv ./ACBDRBXDRB1DQB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop A~C~B DRBX~DRB1~DQB1~DPB1 greedy 0 ACBDRBXDRB1DQB1DPB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./ACBDRBXDRB1DQB1DPB1/pull.$dataset.$pop.txt ./ACBDRBXDRB1DQB1DPB1/glid.$dataset.$pop.txt ./ACBDRBXDRB1DQB1DPB1/freqs.$pop.csv ./ACBDRBXDRB1DQB1DPB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";
  print SLURMIMPUTE "python3 make_blocks_dataset.py $dataset $pop A~C~B DRBX~DRB1~DQB1 greedy 0 ACBDRBXDRB1DQB1 0.01\n";
  print SLURMIMPUTE "python3 impute_GL_dataset.py  $dataset ./ACBDRBXDRB1DQB1/pull.$dataset.$pop.txt ./ACBDRBXDRB1DQB1/glid.$dataset.$pop.txt ./ACBDRBXDRB1DQB1/freqs.$pop.csv ./ACBDRBXDRB1DQB1/impute.$dataset.$pop.csv.gz 0.99 1000\n";


  print SLURMIMPUTE "ELAPSED=\"Elapsed $pop: \$((\$SECONDS / 3600))hrs \$(((\$SECONDS / 60) % 60))min \$((\$SECONDS % 60))sec\"\n";
  print SLURMIMPUTE "echo \$ELAPSED # $pop \n\n\n";

  close SLURMIMPUTE;

}

close SBATCH;