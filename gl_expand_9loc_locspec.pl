#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:  gl_expand_9loc_locspec.pl
# DESCRIPTION:  Reformats 2020 GLID files to expand '/' to '+'
# DATE WRITTEN: March 22, 2020
# WRITTEN BY:   Loren Gragert
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ----------    --------------  -------------------------------------
#
##############################################################################
use strict; # always
use warnings; # always

# my $GLID_file = shift @ARGV;

my $s_loc = shift @ARGV; # selected locus

# load genotype list IDs and fully enumerate genotype lists
my %GLID; # genotype lists by ID
my %GLID_locus; # locus for genotype list
my %glstrings; # unique GL Strings
my %numglids; # number of glids at each locus
my $numgenos = 0; # number of genotypes

# open (GLID,$GLID_file) || die "Missing GLID file $GLID_file\n";
# while (<GLID>) {
while (<STDIN>) {
  chomp;

  my ($glid,$glstring) = split /\,/,$_;

  # print STDERR "GLID: $glid GLString: $glstring\n";
  # store locus for GLID
  if ($glstring =~ m/^\+/) { $glstring = substr($glstring,1) } # some strings start with "+"
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

  if (!defined $GLID_locus{$glid}) {
    print STDERR "GLID with undefined locus $glid\n";
  }

  # skip locus not selected
  if ($GLID_locus{$glid} ne $s_loc) {
    next;
  }

  my @genos_full; # array of genotypes with '+'
  my %unique_genos; # list of unique genotypes

  # GLID #0 is untyped
  if ($glid == 0) { 
    push @genos_full, "UUUU+UUUU";
    print ("0,UUUU+UUUU\n");
    print STDERR "Blank GL String GLID 0 $glstring\n";
    next;
  }

  # make intermediate genotype list
  my @genos_inter = split /\|/, $glstring;
  # make fully enumerated genotype
  foreach my $geno (@genos_inter) {
    my ($l1,$l2) = split /\+/, $geno;
    if (!defined $l2) { $l2 = $l1; }
    if ($l1 eq "") { $l1 = $l2; } # 73354:+DRB1*15:01
    my @l1 = split /\//,$l1;
    my @l2 = split /\//,$l2;

    foreach my $loctyp1 (@l1) {
      foreach my $loctyp2 (@l2) {

        $numgenos++;
        # sort alleles to avoid duplicate genotypes
        # if ($loctyp1 gt $loctyp2) {
        #   ($loctyp1, $loctyp2) = ($loctyp2, $loctyp1);
        # }
        push @genos_full, "$loctyp1+$loctyp2";

        # no dupe genos were found
        # $geno = "$loctyp1+$loctyp2"
        # if (!exists $unique_genos{$geno}) {
        #   $unique_genos{$geno};
        #   push @genos_full, "$loctyp1+$loctyp2";
        # }
        # else {
        #   print STDERR "Dupe geno $geno\n";
        # }
      }
    }
  }

  $glstring = join "|",@genos_full;

  # store locus for GLID
  # if ($glstring =~ m/^A/) { $GLID_locus{$glid} = "A"; }
  # if ($glstring =~ m/^C/) { $GLID_locus{$glid} = "C"; }
  # if ($glstring =~ m/^B/) { $GLID_locus{$glid} = "B"; }
  # if ($glstring =~ m/^DRB1/) { $GLID_locus{$glid} = "DRB1"; }
  # if ($glstring =~ m/^DQB1/) { $GLID_locus{$glid} = "DQB1"; }
  # if ($glstring =~ m/^DRB3/) { $GLID_locus{$glid} = "DRBX"; }
  # if ($glstring =~ m/^DRB4/) { $GLID_locus{$glid} = "DRBX"; }
  # if ($glstring =~ m/^DRB5/) { $GLID_locus{$glid} = "DRBX"; }
  # if ($glstring =~ m/^DRBX/) { $GLID_locus{$glid} = "DRBX"; }
  # if ($glstring =~ m/^DQA1/) { $GLID_locus{$glid} = "DQA1"; }
  # if ($glstring =~ m/^DPB1/) { $GLID_locus{$glid} = "DPB1"; }
  # if ($glstring =~ m/^DPA1/) { $GLID_locus{$glid} = "DPA1"; }

  # $numglids{$GLID_locus{$glid}}++;
  # if (!exists $GLID_locus{$glid}) { print STDERR "GLID has no locus: $_\n"; }

  print "$glid,$glstring\n";

}

# close GLID;
close STDIN;

print STDERR "Number of genotypes in GLID file: $numgenos\n";

exit(0);