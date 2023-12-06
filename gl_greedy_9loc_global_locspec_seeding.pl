#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:  gl_greedy_9loc_locspec.pl
# DESCRIPTION:  Creates a reduced allele list per race and locus starting
#               with common alleles and adding alleles until all donors
#               can be interpreted or until a cutoff is reached.  This results
#               in a reduced genotype list for the EM.
#
# DATE WRITTEN: September 16, 2021
# WRITTEN BY:   Loren Gragert
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ----------    --------------  -------------------------------------
#
##############################################################################
use strict; # always
use warnings; # always

my $s_pop      = shift @ARGV or die "No population provided!\n";

my $s_loc      = shift @ARGV or die "No locus provided!\n";

my $donor_file = shift @ARGV or die "Usage: gl_greedy.pl donor_file glid_file loci hegl_race nemo_race nemo_race_type cutoff cutoff_interp incl_noprobe incl_notyping incl_nohighres";

my $output_dir = shift @ARGV or die "No output directory provided!\n";

my $cfg_dir    = shift @ARGV or die "No cfg directory provided!\n";

my $common_alleles_pop = shift @ARGV or die "No common alleles population provided";

my $cutoff_interp = shift @ARGV or die "No cutoff!\n"; # cutoff for reinterpretation of GLIDs - including UUUU

my $cutoff_interp_classIIA = shift @ARGV or die "No cutoff!\n"; # cutoff for reinterp of DQA1 and DPA1 - 99.99%?

# set to 1 to use initial seed allele list for missing typings
my $cutoff = 1.0; # 1.0 includes all donors - no reason to exclude, just halt reinterpretation
my $common_allele_seeding = 1;
my $missing_typing_initial_list = 0; # for missing typing use seed list and output full list to antigens file
my $missing_typing_cutoff_list = 1; # for missing typing use cutoff list and output that list to antigens file

$| = 1;

my $loci = "A-C-B-DRBX-DRB1-DQA1-DQB1-DPA1-DPB1";
# loci to get from each list
my @locilist = split /-/,$loci;

# load donor glids for each locus
my %ID; # donor typing as set of genotype lists
my %GLID_count; # number of donors with GLID
my $ndonor = scalar %ID;
&loadDonors(\%ID,\%GLID_count,\$ndonor,$donor_file);

print STDERR "$s_pop donors: $ndonor\n";

# load genotype list IDs and fully enumerate genotype lists
my %GLID; # genotype lists by ID
my %GLID_locus; # locus for genotype list
my %numglids; # number of glids at each locus
while (<STDIN>) {
  chomp;

  if ($_ eq "0,UUUU+UUUU") { next;} # skip blank GLID
  if ($_ eq "0,+") { next;} # skip blank GLID
  my ($glid,$glstring) = split /\,/;

  # skip GLIDs that don't exist in cohort
  #if (!exists $GLID_count{$glid}) { next; }

  # store locus for GLID
  if ($glstring =~ m/^A/) { $GLID_locus{$glid} = "A"; }
  if ($glstring =~ m/^C/) { $GLID_locus{$glid} = "C"; }
  if ($glstring =~ m/^B/) { $GLID_locus{$glid} = "B"; }
  if ($glstring =~ m/^DRB1/) { $GLID_locus{$glid} = "DRB1"; }
  if ($glstring =~ m/^DQB1/) { $GLID_locus{$glid} = "DQB1"; }
  if ($glstring =~ m/^DRB3/) { $GLID_locus{$glid} = "DRBX"; }
  if ($glstring =~ m/^DRB4/) { $GLID_locus{$glid} = "DRBX"; }
  if ($glstring =~ m/^DRB5/) { $GLID_locus{$glid} = "DRBX"; }
  if ($glstring =~ m/^DRBX/) { $GLID_locus{$glid} = "DRBX"; }
  if ($glstring =~ m/^DQA1/) { $GLID_locus{$glid} = "DQA1"; }
  if ($glstring =~ m/^DPB1/) { $GLID_locus{$glid} = "DPB1"; }
  if ($glstring =~ m/^DPA1/) { $GLID_locus{$glid} = "DPA1"; }

  # skip locus not selected
  if ($GLID_locus{$glid} ne $s_loc) {
    next;
  }

  $numglids{$GLID_locus{$glid}}++;

  # genotype list
  my @genos = split /\|/, $glstring;

  $GLID{$glid} = [ @genos ];
}


# explicitly handle blank GLID 0
my $glstring;
$glstring = "C*UUUU+C*UUUU";
$GLID_locus{"0"} = "C";

# genotype list
my @genos = split /\|/, $glstring;

$GLID{"0"} = [ @genos ];
$numglids{$GLID_locus{"0"}}++;

$glstring = "DQB1*UUUU+DQB1*UUUU";
$GLID_locus{"999999999"} = "DQB1";
$GLID{"999999999"} = [ @genos ];
$numglids{$GLID_locus{"999999999"}}++;

$glstring = "DRBX*UUUU+DRBX*UUUU";
$GLID_locus{"888888888"} = "DRBX";
$GLID{"888888888"} = [ @genos ];
$numglids{$GLID_locus{"888888888"}}++;

$glstring = "DPA1*UUUU+DPA1*UUUU";
$GLID_locus{"11111111111"} = "DPA1";
$GLID{"11111111111"} = [ @genos ];
$numglids{$GLID_locus{"11111111111"}}++;

$glstring = "DQA1*UUUU+DQA1*UUUU";
$GLID_locus{"55555555"} = "DQA1";
$GLID{"55555555"} = [ @genos ];
$numglids{$GLID_locus{"55555555"}}++;

$glstring = "DPB1*UUUU+DPB1*UUUU";
$GLID_locus{"777777777"} = "DPB1";
$GLID{"777777777"} = [ @genos ];
$numglids{$GLID_locus{"777777777"}}++;

foreach my $locus (@locilist) {

  # skip loci not selected
  if ($locus ne $s_loc) {
    next;
  }

  my $numglids = scalar keys %GLID;
  print STDERR "GLIDs Loaded for $locus: $numglids{$locus}\n";
}



my $finished = 0;
my $cycle = 1;

# load common allele lists or updated allele list after 1st runthrough
my %ALL; # allele list

# load common alleles
&loadAlleleLists(\%ALL,$s_pop,$common_allele_seeding);


# store copy of initial allele list to use for UUUU
my %ALL_init = %ALL;
my %ALL_cutoff = %ALL; # alleles added before hitting cutoff

my %allele_sets; # groups of alleles that were added at the same time
my %stop_reinterp; # per-locus flag for stopping GLID reinterpretation
my %GLID_hits; # number of alleles required to have valid genotype list
my %GLID_reduced; # reduced genotype list
while (!$finished && ($cycle <= 50000)) {

  print STDERR "Cycle: $cycle\n";
  $cycle++;

  # make reduced genotype list
  &makeReducedGL(\%GLID,\%ALL,\%GLID_locus,\%GLID_hits,\%GLID_reduced, \%stop_reinterp);

  # find best single hit allele per locus and update allele lists
  my %doublehits; # doublehits to exclude from ID list
  my %uninterpreted_donors; # donors who have GLID that do not interpret
  &addAlleleSingleHit(\%GLID,\%GLID_count,\%GLID_locus,\%ID,\%ALL,\%ALL_cutoff,\%GLID_hits,\%uninterpreted_donors,\%doublehits,\%allele_sets,$ndonor,$loci,$cutoff,$cutoff_interp,$cutoff_interp_classIIA,\%stop_reinterp,\$finished);

  # output final allele list and genotype list for NEMO
  if ($finished) {
    if ($missing_typing_cutoff_list) { 
      &loadBlankGLstring(\%GLID_reduced,\%ALL_cutoff);
      &outputGL(\%GLID,\%GLID_reduced,\%ALL_cutoff,\%ID,\%uninterpreted_donors,\%doublehits,\%allele_sets,$s_pop);
    }
    elsif ($missing_typing_initial_list) {
      &loadBlankGLstring(\%GLID_reduced,\%ALL_init);
      &outputGL(\%GLID,\%GLID_reduced,\%ALL,\%ID,\%uninterpreted_donors,\%doublehits,\%allele_sets,$s_pop);
    }
    else {
      &loadBlankGLstring(\%GLID_reduced,\%ALL);
      &outputGL(\%GLID,\%GLID_reduced,\%ALL,\%ID,\%uninterpreted_donors,\%doublehits,\%allele_sets,$s_pop);
    }
  }

} # end while

exit (0);



##############################################################################
# Function: loadDonors
##############################################################################
sub loadDonors {

  my ($rID,$rGLID_count,$rndonor,$input_file) = @_;

  print STDERR "input: ",$input_file,"\n";

  my $ndonor_abdrb1_zero = 0;
  # my $ndonor_abdrb1_one = 0;

  open (my $fh_input,"<", $input_file) or die "CANT OPEN FILE $! $0";
  while (<$fh_input>) {
    chomp;
    $_ =~ s/\r//g; # remove DOS carriage return

    # my ($id,$reg_dte,$ctr_cde,$intl_dom_cde,$cntry,$race,$ethnicity,$glid_a,$glid_b,$glid_c,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/,/,$_;
    # my ($id,$race,$ethnicity,$reg_dte,$cntry,$glid_a,$glid_b,$glid_c,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/,/,$_;

    # other pull files have a different format
    my ($id,$race,$ethnicity,$cntry,$glid_a,$glid_b,$glid_c,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1);
    my ($glidlist);
    my ($reg_dte);
    my ($emflag);
    my ($status);
    # print ("Population $s_pop\n");
    # SRTR extracted at Tulane
    if ($s_pop =~ m/^srtr*/) {
      ($id,$glidlist) = split/,/,$_;
      ($glid_a,$glid_c,$glid_b,$glid_drbx,$glid_drb1,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/:/,$glidlist;
    }
    elsif ($s_pop =~ m/^nyukidneyheart*/) {
      ($id,$glidlist) = split/,/,$_;
      ($glid_a,$glid_c,$glid_b,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/:/,$glidlist;
    }
    elsif ($s_pop =~ m/^gustaveroussy*/) {
      ($id,$glidlist) = split/,/,$_;
      ($glid_a,$glid_c,$glid_b,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/:/,$glidlist;
    }
    elsif ($s_pop =~ m/^nyusrtr*/) {
      ($id,$glidlist) = split/,/,$_;
      ($glid_a,$glid_c,$glid_b,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/:/,$glidlist;
    }
    else {
      # NEMO format 2020
      # ($id,$race,$ethnicity,$cntry,$glid_a,$glid_b,$glid_c,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/,/,$_;
      # NEMO format 2022
      ($id,$race,$ethnicity,$reg_dte,$cntry,$emflag,$status,$glid_a,$glid_b,$glid_c,$glid_drb1,$glid_drbx,$glid_dqa1,$glid_dqb1,$glid_dpa1,$glid_dpb1) = split/,/,$_;
    }

    # minimum A B DRB1
    if ($glid_a == 0 || $glid_b == 0 || $glid_drb1 == 0) {
      print STDERR "Missing A B DRB1 GLID 0 $id\n";
      $ndonor_abdrb1_zero++;
      next;
    }
    
    if (($glid_a eq " ") || ($glid_c eq " ") || ($glid_b eq " ") || ($glid_drbx eq " ") || ($glid_drb1 eq " ") || ($glid_dqb1 eq " ")) {
      print STDERR "Donor is missing typing $id $glid_a $glid_c $glid_b $glid_drbx $glid_drb1 $glid_dqb1\n";
      next; 
    }

    # rename dqb1 to "999999999" to avoid mixup with C
    # if ($glid_c == 0) { $glid_c = 1; }
    if ($glid_dqb1 == 0) { $glid_dqb1 = 999999999; }
    if ($glid_drbx == 0) { $glid_drbx = 888888888; } # handle DRBX blanks
    if ($glid_dqa1 == 0) { $glid_dqa1 = 55555555; } 
    if ($glid_dpa1 == 0) { $glid_dpa1 = 11111111111; } 
    if ($glid_dpb1 == 0) { $glid_dpb1 = 777777777; } 

    if ($glid_dqb1 == 1) { $glid_dqb1 = 999999999; }
    if ($glid_drbx == 1) { $glid_drbx = 888888888; } # handle DRBX blanks
    if ($glid_dqa1 == 1) { $glid_dqa1 = 55555555; } 
    if ($glid_dpa1 == 1) { $glid_dpa1 = 11111111111; } 
    if ($glid_dpb1 == 1) { $glid_dpb1 = 777777777; } 

    my $glid_locus = "";
    if ($s_loc eq "A") {
      $glid_locus = $glid_a;
    }
    elsif ($s_loc eq "C") {
      $glid_locus = $glid_c;
    }
    elsif ($s_loc eq "B") {
      $glid_locus = $glid_b;
    }
    elsif ($s_loc eq "DRBX") {
      $glid_locus = $glid_drbx;
    }
    elsif ($s_loc eq "DRB1") {
      $glid_locus = $glid_drb1;
    }
    elsif ($s_loc eq "DQA1") {
      $glid_locus = $glid_dqa1;
    }
    elsif ($s_loc eq "DQB1") {
      $glid_locus = $glid_dqb1;
    }
    elsif ($s_loc eq "DPA1") {
      $glid_locus = $glid_dpa1;
    }
    elsif ($s_loc eq "DPB1") {
      $glid_locus = $glid_dpb1;
    }
    else {
    }

    # switch to genomic order
    $$rID{$id} = $glid_locus;
    $$rGLID_count{$glid_locus}++;
  

    $$rndonor++;
  }
  close $fh_input;
  
print STDERR "Missing A or B or DRB1 GLID 0 $ndonor_abdrb1_zero\n";
# print STDERR "Missing A or B or DRB1 GLID 1 $ndonor_abdrb1_one\n";

} # end sub loadDonors


##########################################################################
# Function: loadAlleleLists - loads allele lists - Initial and Updated
##########################################################################
sub loadAlleleLists {

  my ($rALL, $nemo_race, $common_allele_seeding) = @_;
  
  # my %h_classI = (
  #   "A" => 1,
  #   "C" => 1,
  #  "B" => 1
  # );

  # reference alleles from IMGT/HLA for generating first pass 
  $$rALL{"A*01:01"}++;
  $$rALL{"B*07:02"}++;
  $$rALL{"C*01:02"}++;
  $$rALL{"DRB1*01:01"}++;
  if ($nemo_race eq "NG") { 
    $$rALL{"DRB1*03:01"}++; # needed for NG WMDA population to interpret genotype
  }
  $$rALL{"DRB3*01:01"}++;
  $$rALL{"DRB4*01:01"}++;
  $$rALL{"DRB5*01:01"}++;
  $$rALL{"DQA1*01:01"}++;
  $$rALL{"DQB1*05:01"}++;
  $$rALL{"DPA1*01:03"}++;
  $$rALL{"DPB1*01:01"}++;

  foreach my $locus (@locilist) {


    # handle blank allele
    $$rALL{"$locus*UUUU"}++;

    if ($common_allele_seeding) {

      my $allele_filename = "";
      $allele_filename = $cfg_dir . "/common_alleles_" . $common_alleles_pop ."_" . $locus . ".cfg";
      # if (($s_pop eq "RS") || ($s_pop eq "HU")) {
      #   $allele_filename = $cfg_dir."/common_alleles_CAU_" . $locus . ".cfg";
      # }
      # if ($s_pop eq "GLOBAL") {
      #   $allele_filename = $cfg_dir."/common_alleles_GLOBAL_" . $locus . ".cfg";
      # }
      # if ($s_pop eq "ukt") {
      #   $allele_filename = $cfg_dir."/common_alleles_CAU_" . $locus . ".cfg";
      # }
      # if ($s_pop eq "ukt") {
      #   $allele_filename = $cfg_dir."/common_alleles_CAU_" . $locus . ".cfg";
      # }

      open (ALLELE_LIST,"<", "$allele_filename") or die "$allele_filename Cant open $!\n";

      my $numalleles = 0;
      while (<ALLELE_LIST>) {
        chomp;
        my ($allele,$count,$freq) = split /,/,$_;
        if ($allele eq "Allele") { next; }
        $$rALL{$allele}++;
        $numalleles++;
      }

      $numalleles++;

      close ALLELE_LIST;
      
      print STDERR "Loaded $numalleles $locus alleles from $allele_filename\n";

    } # end if common_allele_seeding

  } # end foreach locus

} # end sub loadAlleleLists

##########################################################################
# Function: makeReducedGL - Update enumerated GLIDs to reduced allele list
##########################################################################
sub makeReducedGL {

  my ($rGLID,$rALL,$rGLID_locus,$rGLID_hits,$rGLID_reduced,$rstop_reinterp) = @_;

  my $genos_orig = 0; # total number of genotypes for all GLIDs 
  my $genos_reduced = 0; # number of genotypes after reduction to allele list
  foreach my $glid (keys %$rGLID) {
    my $locus = $$rGLID_locus{$glid}; 

    # skip GLIDs not from selected locus
    if ($locus ne $s_loc) {
      next;
    }

    # skip GLID if stop_reinterp flag has been set and GLID has list
    if ($$rstop_reinterp{$locus}) {
      if (exists $$rGLID_hits{$glid}) {
        if ($$rGLID_hits{$glid} == 0) {
         # print STDERR "Skipped reinterp of $glid\n";
          next;
        }
      }
    }

    my @newgenos; # genotypes after reduced allele list

    # fewest number of alleles required for genotype list entry
    my $min_hits = 2; 

    for my $i (0 .. $#{ $$rGLID{$glid} } ) {
      
	    my $geno = $$rGLID{$glid}[$i];
	  
      $genos_orig++;

      my ($loctyp1,$loctyp2) = split /\+/, $geno;
	    # print STDERR "$glid $geno\n" if !defined $loctyp1 || $loctyp1 !~ /\S/ || !defined $loctyp2 || $loctyp2 !~ /\S/;
	  
      my ($loc1,$typ1) = split /\*/,$loctyp1;
      my ($loc2,$typ2) = split /\*/,$loctyp2;
      

      my $hits = 0; # hits are the number of alleles required to make a geno
      if (!exists $$rALL{"$loc1*$typ1"}) {
        $hits++;
      }
      if (!exists $$rALL{"$loc2*$typ2"}) {
        $hits++;
      }
      
      # print "$glid $loctyp1 $loctyp2 $hits hits\n";
      
      if (exists $$rALL{"$loc1*$typ1"} && exists $$rALL{"$loc2*$typ2"}) {
        # sort genotype so lowest number allele is first
        if ($loctyp1 gt $loctyp2) {
          my $tmp = $loctyp1;
          $loctyp1 = $loctyp2;
          $loctyp2 = $tmp;
          $geno = "$loctyp1+$loctyp2";
        }
        push(@newgenos,$geno);  # both alleles exist - add to new genotype list
        $genos_reduced++;
      }

      if ($hits < $min_hits) { $min_hits = $hits; }

    }

    # store reduced genotype list
    $$rGLID_reduced{$glid} = [ @newgenos ];

    # record how many alleles this typing needs to have a valid genotype
    $$rGLID_hits{$glid} = $min_hits;

  } # end foreach glid

  my $reduction_factor = $genos_orig/$genos_reduced;
  print STDERR "Genotype List Reduction: $genos_orig to $genos_reduced - Factor of $reduction_factor\n";

} # end sub makeReducedGL


##########################################################################
# Function: addAlleleSingleHit - add allele to list based on single hits
##########################################################################
sub addAlleleSingleHit {

  my ($rGLID,$rGLID_count,$rGLID_locus,$rID,$rALL,$rALL_cutoff,$rGLID_hits,$runinterpreted_donors,$rdoublehits,$rallele_sets,$rndonor,$loci,$rcutoff,$rcutoff_interp,$rcutoff_interp_classIIA,$rstop_reinterp,$rfinished) = @_;


  my %ALL_hits; # count number of donors with hits per allele
  my %ALL_glids; # list of GLIDs where allele would give interpretation
  my %ndonor_valid; # number of donors with glid
  my %numglids_valid; # number of glids interpreted
  my $ndonor_doublehit = 0; # number of donors with double hit glids
  foreach my $glid (keys %$rGLID) {

    # handle case where no donors from population have a certain GLID
    if (!exists $$rGLID_count{$glid}) {$$rGLID_count{$glid} = 0;}
    
    my $ndonor_GLID = $$rGLID_count{$glid}; # number of donors with GLID
    my $locus = $$rGLID_locus{$glid};

    # skip locus that was not selected
    if ($locus ne $s_loc) {
      next;
    }

    # print STDERR "glid: $glid loc: $locus count: $$rGLID_count{$glid} hits: $$rGLID_hits{$glid}\n";

    # count valid glids and donors
    if ($$rGLID_hits{$glid} == 0) { 
      $ndonor_valid{$locus} += $ndonor_GLID;
      $numglids_valid{$locus}++;
    }

    # count number of donors with double hit GLIDs
    if ($$rGLID_hits{$glid} == 2) {
      $ndonor_doublehit += $$rGLID_count{$glid};
      if ($$rGLID_count{$glid} > 0) {
        $$rdoublehits{$glid}++;
      }
    }

    # single-hit GLIDs only
    if ($$rGLID_hits{$glid} != 1) { next; }

    my %geno_counted; # flag for if genotype was already counted for glid
    # count number of donors who need single allele to have valid list
    for my $i (0 .. $#{ $$rGLID{$glid} } ) {
      my $geno = $$rGLID{$glid}[$i];
      # print STDERR "$glid $geno\n";
      my ($loctyp1,$loctyp2) = split /\+/, $geno;
      my ($loc1,$typ1) = split /\*/,$loctyp1;
      my ($loc2,$typ2) = split /\*/,$loctyp2;
     
      # skip double hits
      if (!exists $$rALL{"$loc1*$typ1"} && !exists $$rALL{"$loc2*$typ2"}) {
        # print STDERR "Double Hit $glid $loctyp1 $loctyp2\n";
         next;
      }
          
      # count single hits
      if (!exists $$rALL{"$loc1*$typ1"}) {
       # print STDERR "Single Hit $glid $loctyp1\n";
        if (!exists $geno_counted{$geno}) {
          $ALL_hits{"$loc1*$typ1"} += $ndonor_GLID;
          push @{$ALL_glids{"$loc1*$typ1"}}, $glid;
          $geno_counted{$geno}++;
        }
      }
      if (!exists $$rALL{"$loc2*$typ2"}) {
        if (!exists $geno_counted{$geno}) {
          $ALL_hits{"$loc2*$typ2"} += $ndonor_GLID;
          push @{$ALL_glids{"$loc2*$typ2"}}, $glid;
          $geno_counted{$geno}++;
        }
        #print STDERR "Single Hit $glid $loctyp2\n";
      }
    }

  } # end foreach glid

  # set flag to stop reinterp
  foreach my $locus (@locilist) {

    if ($locus ne $s_loc) {
      next;
    }

    my $percent_donors_locus = $ndonor_valid{$locus} / $rndonor;
    # print STDERR "$locus $rndonor $ndonor_valid{$locus}\n";
    print STDERR "$ndonor_valid{$locus} of $rndonor donors ($percent_donors_locus) at $locus interpreted\n";
    print STDERR "$numglids_valid{$locus} of $numglids{$locus} GLIDs at $locus interpreted\n";
    # set stop_reinterp flag if enough donors have interpreted
    if (($locus eq "DQA1" || $locus eq "DPA1")) {
      if ($percent_donors_locus >= $rcutoff_interp_classIIA) {
        print STDERR "GLID Reinterpretation halted at cutoff $rcutoff_interp_classIIA for locus $locus\n";
        $$rstop_reinterp{$locus} = 1;
      }
    }
    else {
      if ($percent_donors_locus >= $rcutoff_interp) {
        print STDERR "GLID Reinterpretation halted at cutoff $rcutoff_interp for locus $locus\n";
        $$rstop_reinterp{$locus} = 1;
      }
    }
  }

  print STDERR "$ndonor_doublehit donors with double hit GLIDs\n";
  # foreach my $glid (keys %$rdoublehits) {
    # print "Donor has Double Hit $glid\n";
  # }


  my $total_valid_donors = 0; # donors with GLIDs at all loci interpreting
  my $ID_key_count = scalar keys %$rID;
  # print STDERR "Unique IDs: $ID_key_count\n";
  foreach my $id (keys %$rID) {
    my $glid_locus = $ID{$id};

    if (!exists($$rGLID_hits{$glid_locus})) {
      print STDERR "$s_loc GLID not defined: ID $id GLID $glid_locus\n";
    }

    if ($$rGLID_hits{$glid_locus} == 0) {
      $total_valid_donors++;
    }
    else {
      $$runinterpreted_donors{$id}++;
    }

    #print STDERR "Glid_hits -> ",join(":",$id,$$rGLID_hits{$glid_a},$$rGLID_hits{$glid_c},$$rGLID_hits{$glid_b},
    #  $$rGLID_hits{$glid_drbx},$$rGLID_hits{$glid_drb1},$$rGLID_hits{$glid_dqb1},
    #  $$rGLID_hits{$glid_dpa1},$$rGLID_hits{$glid_dqa1},$$rGLID_hits{$glid_dpb1}),"\n";

  }

  my $num_invalid_donors = scalar keys %$runinterpreted_donors;
  if ($num_invalid_donors <= 10) {
    foreach my $id (keys %$runinterpreted_donors) {
      print STDERR "$id not interpreted $ID{$id}\n";
    }
  }
  
  my $percent_valid_donors = $total_valid_donors / $rndonor;
  print STDERR "$total_valid_donors of $rndonor donors ($percent_valid_donors) interpreted at all loci\n";


  # set finished flag if enough donors interpret
  if (($percent_valid_donors >= $rcutoff) || (($total_valid_donors + $ndonor_doublehit) == $rndonor)) {
    print STDERR "Met Cutoff - Finished\n";
    $$rfinished = 1;
    return;
  }

  my %loc_added; # flags for if allele for that locus was found
  foreach my $locus (@locilist) {
    my $prev_numhits = 0; # number of hits in last round
    my $prev_typ = "UUUU"; # last allele checked
    my $prev_loc = "HLA";
    my %allele_dupe_hits; # store all alleles with same number of hits

    # skip locus that was not selected
    if ($locus ne $s_loc) {
      next;
    }

    foreach my $loctyp (sort {$ALL_hits{$b} <=> $ALL_hits{$a}} keys %ALL_hits) {
      my $numhits = $ALL_hits{$loctyp};

      my ($loc, $typ) = split /\*/, $loctyp;
  
      if (($loc eq "DRB3") || ($loc eq "DRB4") || ($loc eq "DRB5")) {
        if ($locus ne "DRBX") { next; }
      }
      else {
        if ($loc ne $locus) { next; }
      }

      # stop adding at locus when there are no more alleles
      if ($numhits == 0) { next; }

      # skip alleles for loci already seen
      if (exists $loc_added{$loc}) { next; } 

      # print STDERR "$loctyp has $numhits hits\n";
    
      if ($numhits == $prev_numhits) { # check for duplicate hits
        # print STDERR "Duplicate hits $numhits $prev_numhits $loctyp\n";
        $prev_numhits = $numhits;
        $prev_typ = $typ;
        $prev_loc = $loc;
        $allele_dupe_hits{$loctyp} = 1;
      }
      elsif ($prev_typ eq "UUUU") {  # do nothing on first pass
        # print STDERR "UUUU initial $loctyp\n";
        $prev_numhits = $numhits;
        $prev_typ = $typ;
        $prev_loc = $loc;
        $allele_dupe_hits{$loctyp} = 1;
      }      
      else { # skip all cases when numbers of hits decrements
        # print STDERR "Hits decremented $loctyp\n";
        $loc_added{$locus}++;
        next;
      }
    }

    my $dupe_numhits = scalar keys %allele_dupe_hits;

    # add last allele if no other allele has the same number of hits
    if ($dupe_numhits == 1) {
      print STDERR "$prev_loc*$prev_typ - $prev_numhits hits added\n";
      $$rALL{"$prev_loc*$prev_typ"}++; # add allele to list
      if (!$$rstop_reinterp{$locus}) {
        $$rALL_cutoff{"$prev_loc*$prev_typ"}++; # add allele to list before cutoff
      }
      $loc_added{$locus}++; # flag for locus
      next;
    }
    else {
      if ($prev_numhits != 0) {
        print STDERR "Number of $locus alleles with same number of hits: $dupe_numhits\n";
      }
    }

    # compare sorted list of uninterpreted GLIDs where alleles appear
    my %allele_glid_matches; # list of alleles with same GLIDs
    my %allele_glid_nummatches; # list of alleles with same GLIDs
    foreach my $loctyp1 (sort keys %allele_dupe_hits) {
      my @glidlist1 = @{$ALL_glids{$loctyp1}};
      foreach my $loctyp2 (sort keys %allele_dupe_hits) {
        my @glidlist2 = @{$ALL_glids{$loctyp2}};
        # compare arrays
        my $are_equal = compare_arrays(\@glidlist1, \@glidlist2);
        if ($are_equal) {
          push @{$allele_glid_matches{$loctyp1}},$loctyp2;
          $allele_glid_nummatches{$loctyp1}++;
          # print STDERR "$loctyp1 and $loctyp2 GLID arrays are equal\n";
        }
      }
    }

    # add all alleles where GLID lists are the same for largest group with
    # matching GLID lists
    my $first_allele_group = 1;
    my @allele_set; # set of alleles added at the same time
    foreach my $loctyp1 (sort {$allele_glid_nummatches{$b} <=> $allele_glid_nummatches{$a}} keys %allele_glid_nummatches) {
      # only do first group of alleles
      if ($first_allele_group == 0) { next; } 
      # add all alleles that travel together in same group of GLIDs
      foreach my $loctyp2 (@{$allele_glid_matches{$loctyp1}}) {
        $$rALL{$loctyp2}++; # add allele to list
        if (!$$rstop_reinterp{$locus}) {
          $$rALL_cutoff{"$loctyp2"}++; # add allele to list before cutoff
        }
        push @allele_set, $loctyp2;
        print STDERR "$loctyp2 - $prev_numhits hits added\n";
      }
      $first_allele_group = 0;
      if (scalar @allele_set > 1) {
        my $allele_group = join "+", @allele_set;
        $$rallele_sets{$allele_group} = $prev_numhits;
        # print STDERR "Group $allele_group recorded\n";
      }
    }

  } # end locus 

  if (0) {
  # version that adds lowest numbered allele first
  # add most prevalent allele per locus
  my %loc_added; # flags for if allele for that locus was found
  foreach my $locus (@locilist) {
    my $prev_numhits = 0; # number of hits in last round
    my $prev_typ = "UUUU"; # last allele checked
    foreach my $loctyp (sort {$ALL_hits{$b} <=> $ALL_hits{$a}} keys %ALL_hits) {
      my $numhits = $ALL_hits{$loctyp};

      my ($loc, $typ) = split /\*/, $loctyp;
  
      if ($loc ne $locus) { next; }

      # stop adding at locus when there are no more alleles
      if ($numhits == 0) { next; }

      # skip alleles for loci already seen
      if (exists $loc_added{$loc}) { next; } 

      # print STDERR "$loctyp has $numhits hits\n";
    
      if ($numhits == $prev_numhits) { # check for duplicate hits
        # next time, load lowest numbered allele
        if ($prev_typ > $typ) {
          $prev_typ = $typ;
          $prev_numhits = $numhits;
        }      
      }
      elsif ($prev_typ eq "UUUU") {  # do nothing on first pass
        $prev_typ = $typ;
        $prev_numhits = $numhits;
      }
      # current allele has more hits - add lowest numbered allele seen so far
      else {
        $$rALL{"$loc*$prev_typ"}++; # add allele to list
        $loc_added{$loc}++; # flag for locus

        print STDERR "$loc*$prev_typ - $prev_numhits hits added\n";
      }
      

    } # end loctyp

    # handle one hit case
    if ($prev_numhits == 1) {
      $$rALL{"$locus*$prev_typ"}++; # add allele to list
      print STDERR "$locus*$prev_typ - $prev_numhits hits added\n";      
    }

  } # end locus 

  } # end if 0 - lowest allele first


 
  
} # end sub addAlleleSingleHit


##########################################################################
# Function: loadBlankGLstring - Load GLstring for blank GLID
##########################################################################
sub loadBlankGLstring {

  my ($rGLID_reduced,$rALL) = @_;

  my @newgenos_C; # genotype list for no typing at locus
  my @newgenos_DQB1;
  my @newgenos_DRBX;
  my @newgenos_DQA1;
  my @newgenos_DPA1;
  my @newgenos_DPB1;

  foreach my $loctyp1 (keys %$rALL) {
    my ($loc1,$typ1) = split /\*/,$loctyp1;
    if ($loc1 ne "C") { next; }
    if ($typ1 eq "UUUU") { next; }
    foreach my $loctyp2 (keys %$rALL) {
      my ($loc2,$typ2) = split /\*/,$loctyp2;
      if ($loc2 ne "C") { next; }
      if ($typ2 eq "UUUU") { next; }
      my $geno = "$loctyp1+$loctyp2";
      push (@newgenos_C,$geno);
      # print "C initial genos $geno $loc1 $loc2\n";
    }
  }
  $$rGLID_reduced{"0"} = [ @newgenos_C ];


  foreach my $loctyp1 (keys %$rALL) {
    my ($loc1,$typ1) = split /\*/,$loctyp1;
    if ($loc1 ne "DQB1") { next; }
    if ($typ1 eq "UUUU") { next; }
    foreach my $loctyp2 (keys %$rALL) {
      my ($loc2,$typ2) = split /\*/,$loctyp2;
      if ($loc2 ne "DQB1") { next; }
      if ($typ2 eq "UUUU") { next; }
      my $geno = "$loctyp1+$loctyp2";
      push (@newgenos_DQB1,$geno);
    }
  }
  $$rGLID_reduced{"999999999"} = [ @newgenos_DQB1 ];

  foreach my $loctyp1 (keys %$rALL) {
    my ($loc1,$typ1) = split /\*/,$loctyp1;
    if ($loc1 ne "DQA1") { next; }
    if ($typ1 eq "UUUU") { next; }
    foreach my $loctyp2 (keys %$rALL) {
      my ($loc2,$typ2) = split /\*/,$loctyp2;
      if ($loc2 ne "DQA1") { next; }
      if ($typ2 eq "UUUU") { next; }
      my $geno = "$loctyp1+$loctyp2";
      push (@newgenos_DQA1,$geno);
    }
  }
  $$rGLID_reduced{"55555555"} = [ @newgenos_DQA1 ];

    foreach my $loctyp1 (keys %$rALL) {
    my ($loc1,$typ1) = split /\*/,$loctyp1;
    if ($loc1 ne "DPB1") { next; }
    if ($typ1 eq "UUUU") { next; }
    foreach my $loctyp2 (keys %$rALL) {
      my ($loc2,$typ2) = split /\*/,$loctyp2;
      if ($loc2 ne "DPB1") { next; }
      if ($typ2 eq "UUUU") { next; }
      my $geno = "$loctyp1+$loctyp2";
      push (@newgenos_DPB1,$geno);
    }
  }
  $$rGLID_reduced{"777777777"} = [ @newgenos_DPB1 ];

  
    foreach my $loctyp1 (keys %$rALL) {
    my ($loc1,$typ1) = split /\*/,$loctyp1;
    if ($loc1 ne "DPA1") { next; }
    if ($typ1 eq "UUUU") { next; }
    foreach my $loctyp2 (keys %$rALL) {
      my ($loc2,$typ2) = split /\*/,$loctyp2;
      if ($loc2 ne "DPA1") { next; }
      if ($typ2 eq "UUUU") { next; }
      my $geno = "$loctyp1+$loctyp2";
      push (@newgenos_DPA1,$geno);
    }
  }
  $$rGLID_reduced{"11111111111"} = [ @newgenos_DPA1 ];

  # add all allele combos to genotype list
  foreach my $loctyp1 (keys %$rALL) {
    my ($loc1,$typ1) = split /\*/,$loctyp1;
    if (($loc1 eq "DRB3") || ($loc1 eq "DRB4") || ($loc1 eq "DRB5")) {
      foreach my $loctyp2 (keys %$rALL) {
        my ($loc2,$typ2) = split /\*/,$loctyp2;
        if (($loc2 eq "DRB3") || ($loc2 eq "DRB4") || ($loc2 eq "DRB5")) {
          my $geno = "$loctyp1+$loctyp2";
          push (@newgenos_DRBX,$geno);
        }
      }
      # add blank option for locus 2
      my $geno = "$loctyp1+DRBX*NNNN";        
      push (@newgenos_DRBX,$geno);
    }
  }

  # add blank option for both loci
  my $geno = "DRBX*NNNN+DRBX*NNNN";        
  push (@newgenos_DRBX,$geno);

  $$rGLID_reduced{"888888888"} = [ @newgenos_DRBX ];

  # add UUUU option to each homozygous genotype
  foreach my $glid (keys %$rGLID_reduced) {
    my $glstring = join "|",@{$$rGLID_reduced{$glid}};
    if (!(($glstring =~ m/^DRB3/) || ($glstring =~ m/^DRB4/) || ($glstring =~ m/^DRBX/)
          || ($glstring =~ m/^DRB5/))) { next; } # skip non-DRBX GLIDs 
    my @geno = @{$$rGLID_reduced{$glid}};
    my $homo_DRBX = 0; # 1 if homozygous typing found
    my %loctyp1; # list of DRBX typings
    foreach my $geno (@geno) {
      my ($loctyp1,$loctyp2) = split /\+/,$geno;
      if ($loctyp1 eq $loctyp2) { $homo_DRBX = 1; }
      $loctyp1{$loctyp1} = 1; 
    }
    if ($homo_DRBX == 1) {
      foreach my $loctyp1 (keys %loctyp1) {
          $geno = "$loctyp1+DRBX*NNNN";
          push (@geno,$geno);
          my ($loc1, $typ1) = split /\*/,$loctyp1;
          if (($loc1 eq "DRB3") || ($loc1 eq "DRB5")) {
            $geno = "$loctyp1+DRB4*01:01";
            push (@geno,$geno);
            $geno = "$loctyp1+DRB4*01:03";
            push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*01:02";
            # push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*01:04";
            # push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*01:05";
            # push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*01:07";
            # push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*01:08";
            # push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*02:01N";
            # push (@geno,$geno);
            # $geno = "$loctyp1,DRB4*03:01N";
            # push (@geno,$geno);
          }
          $$rGLID_reduced{"$glid"} = [ @geno ];
          
        } # end foreach loctyp1
      
    } # end if homo_DRBX
  } # end foreach GLID
 


} # end sub loadBlankGLstring


##########################################################################
# Function: outputGL - output reduced genotype list for new glid extract
#                      and extract of interpretable donors for nemo
##########################################################################
sub outputGL {

  my ($rGLID,$rGLID_reduced,$rALL,$rID,$runinterpreted_donors,$rdoublehits,$rallele_sets,$nemo_race) = @_;

  my $nemo_race_type = "DETAILED_ETH";
  # print reduced GLID file 
  my $glid_filename = $output_dir."/glid.$nemo_race.$s_loc.txt";

  print STDERR "Printing reduced glid file $glid_filename\n";

  open (my $fh_glid,">","$glid_filename");
  # my %GLID_number; # store by glstring to avoid duplicate GLIDs 
  # foreach my $glid (keys %$rGLID_reduced) {
  #   my $glstring = join "|",@{$$rGLID_reduced{$glid}};
  #   $GLID_number{$glstring} = $glid;
  #   # print GLID "$glid:$glstring\n";
  # }
  #my %h = map{ join "|",@{$$rGLID_reduced{$_}} => $_ } (keys %$rGLID_reduced);
  my %GLID_number;#   = map{ my $gl = join "|",@{$$rGLID_reduced{$_}};$gl => $_ } (keys %$rGLID_reduced);
  my %GLID_glstring;# = map{ my $gl = join "|",@{$$rGLID_reduced{$_}};$_  => $gl } (keys %$rGLID_reduced);
  foreach my $glid (keys %$rGLID_reduced){
      my $gl = join "|",@{$$rGLID_reduced{$glid}};
      $GLID_number{$gl}     = $glid;
      $GLID_glstring{$glid} = $gl;
  }

  foreach my $glstring (keys %GLID_number) {
    print $fh_glid "$GLID_number{$glstring},$glstring\n" if defined $glstring && $glstring =~ /\S/;
  }
  close $fh_glid;

  my $num_reduced_glids = scalar keys %$rGLID_reduced;
  my $num_final_glids   = scalar keys %GLID_number;
  print STDERR "Compressing duplicate glstrings $num_reduced_glids to $num_final_glids\n"; 

  # output new donor file with rolled up GLIDs
  my $donor_filename = $output_dir."/pull.$nemo_race.$s_loc.txt";

  print STDERR "Printing reduced donor file $donor_filename\n";

  open (my $fh_donor,">","$donor_filename");
  foreach my $id (keys %$rID) {
    if (exists $$runinterpreted_donors{$id}) { next; }  # skip uninterpreted
    #my $glids = $$rID{$id};
    #my @glids = split /:/, $glids;
    # my @newglids; # GLIDs after duplicate glstring compression
    # foreach(split /:/,$$rID{$id}) {
    #   #my $glstring = join "|",@{$$rGLID_reduced{$glid}};
    #   #my $new_glid = $GLID_number{$glstring};
    #   #push (@newglids,$new_glid);
    #   push (@newglids,$GLID_number{join "|",@{$$rGLID_reduced{$_}}});
    # }
    #my $newglids  = join ":",@newglids;
    my $glid = $$rID{$id};
    my $glstring = $GLID_glstring{$glid};
    my $newglid = $GLID_number{$glstring};
    print $fh_donor "$id,$newglid\n";
  }
  close $fh_donor;

  # output donor file with doublehits
  my $doublehit_filename = $output_dir."/doublehit.$nemo_race.$s_loc.txt";

  print STDERR "Printing list of doublehit GLIDs $doublehit_filename\n";

  open (my $fh_double,">","$doublehit_filename");
  foreach my $glid (keys %$rGLID) {
    if (exists $$rdoublehits{$glid}) {
      my $glstring = join "|",@{$$rGLID{$glid}};
      print $fh_double "$glid,$glstring\n";
    }
  }
  close $fh_double;

  # output donor file with uninterpreted donors
  my $uninterp_filename = $output_dir."/uninterp.$nemo_race.$s_loc.txt";

  print STDERR "Printing list of uninterpreted donors $uninterp_filename\n";

  open (my $fh_uninterp,">","$uninterp_filename");
  foreach my $id (keys %$rID) {
    if (exists $$runinterpreted_donors{$id}) {
      print $fh_uninterp "$id,$$rID{$id}\n";
    }
  }
  close $fh_uninterp;


  # output allele list - either pre-cutoff or full allele list depending on mode
  foreach my $locus (@locilist) {

    # skip loci not selected
    if ($locus ne $s_loc) {
      next;
    }

    my $antigens_filename = $output_dir."antigens.$nemo_race.$s_loc.txt";

    print STDERR "Printing allele list $antigens_filename\n";

    open (my $fh_antigen,">","$antigens_filename");
    foreach my $loctyp (sort keys %$rALL) {
      my ($loc,$typ) = split /\*/,$loctyp;
      if (($loc eq "DRB3") || ($loc eq "DRB4") || ($loc eq "DRB5")) {
        $loc = "DRBX";
      }
      if ($locus ne $loc) { next; }
      if ($typ eq "UUUU") { next; }
      print $fh_antigen "$loctyp\n";
    }
    close $fh_antigen;
  }

  # output list of groups of alleles added at the same hits
  my $allelesets_filename = $output_dir."allelesets.$nemo_race.$s_loc.txt";

  print STDERR "Printing allele sets $allelesets_filename\n";

  open (my $fh_alleleset,">","$allelesets_filename");
  foreach my $allele_list (sort keys %$rallele_sets) {
    # print STDERR "$allele_list\n";
    print $fh_alleleset "$allele_list;$$rallele_sets{$allele_list}\n";
  }
  close $fh_alleleset;


} # end sub outputGL


sub compare_arrays {
	my ($first, $second) = @_;
	no warnings;  # silence spurious -w undef complaints
	return 0 unless @$first == @$second;
	for (my $i = 0; $i < @$first; $i++) {
	    return 0 if $first->[$i] ne $second->[$i];
	}
	return 1;
}  
