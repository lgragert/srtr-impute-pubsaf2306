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


my $sbatch_file = "run_sbatch_greedy_srtr_loni.sh";
open (SBATCH, ">$sbatch_file") || die "Missing $sbatch_file\n";


foreach my $pop (@pops) {

  my $pop_greedy = $pop;

  my $slurm_greedy_file = "./slurm/run_slurm_srtr_" . $pop . "_greedy_loni.sh";
  open (SLURMGREEDY, ">$slurm_greedy_file") || die "Missing $slurm_greedy_file\n";

  print SBATCH "jid_greedy_$pop=\$(sbatch $slurm_greedy_file)\n";
  print SBATCH "echo \$jid_greedy_$pop greedy $pop\n";
  print SBATCH "echo \${jid_greedy_$pop\#\#\* }\n";

  my $job_greedy = $pop . "_greedy";


  print SLURMGREEDY "#!/bin/bash\n";
  print SLURMGREEDY "#SBATCH --job-name=$job_greedy   ### Job Name\n";
  print SLURMGREEDY "#SBATCH --output=./logs/srtr_$pop.greedy.out       ### File in which to store job output\n";
  print SLURMGREEDY "#SBATCH --error=./logs/srtr_$pop.greedy.err        ### File in which to store job error messages\n";
  print SLURMGREEDY "#SBATCH --account=loni_hla_em03   ### Allocation\n";
  print SLURMGREEDY "#SBATCH --partition=single    # Partition default=single\n";
  # print SLURMGREEDY "#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)\n";
  print SLURMGREEDY "#SBATCH --time=0-24:00:00     ### Wall clock time limit in Days-HH:MM:SS\n";
  # print SLURMGREEDY "#SBATCH --mem=128000          ### Memory in MB (default memory is 64000 - 64GB)\n";
  print SLURMGREEDY "#SBATCH --nodes=1             ### Node count required for the job\n";
  print SLURMGREEDY "#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node\n\n\n";

  print SLURMGREEDY "SECONDS=0 # $pop\n\n";

  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl A | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop A ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl B | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop B ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl C | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop C ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl DRBX | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop DRBX ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl DRB1 | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop DRB1 ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl DQA1 | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop DQA1 ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl DQB1 | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop DQB1 ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl DPA1 | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop DPA1 ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "cat ./pull/glid.srtr.$pop.txt | perl gl_expand_9loc_locspec.pl DPB1 | perl gl_greedy_9loc_global_locspec_seeding.pl srtr.$pop DPB1 ./pull/pull.srtr.$pop.txt ./greedy/ ./cfg/ $pop_greedy 0.999 0.9999\n";
  print SLURMGREEDY "python3 combine_greedy_locfiles.py srtr.$pop\n";
  
  print SLURMGREEDY "ELAPSED=\"Elapsed $pop: \$((\$SECONDS / 3600))hrs \$(((\$SECONDS / 60) % 60))min \$((\$SECONDS % 60))sec\"\n";
  print SLURMGREEDY "echo \$ELAPSED # $pop \n\n\n";

  close SLURMGREEDY;

}

close SBATCH;