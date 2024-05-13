import pandas as pd
import re
from aa_matching_msf_genie_MOD import *

eplet_DQ = pd.read_csv("Eplet_Registry_DQ.txt", sep='\t') 

DQ_eplet_positions_all = {}

for eplet_polymorphism in eplet_DQ['Polymorphic']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_all[position] = 1

# Created a dictionary of eplet positions listed in Eplet_Registry_DQ.txt
# Key of AA position
# Value of 1 if it has been seen
DQ_eplet_positions_all

# Count number of mismatches at position between donor and recip for DQ
def count_AA_Mismatches_DQ(self, aa1_donor,aa2_donor,aa3_donor,aa4_donor,aa1_recip,aa2_recip,aa3_recip,aa4_recip):
  mm_count = 0
  if (aa1_donor != aa1_recip):
    mm_count+=1
  if (aa2_donor != aa2_recip):
    mm_count+=1
  if (aa1_donor != aa2_recip):
    mm_count+=1
  if (aa2_donor != aa1_recip):
    mm_count+=1
  if (aa3_donor != aa3_recip):
    mm_count+=1
  if (aa4_donor != aa4_recip):
    mm_count+=1
  if (aa3_donor != aa4_recip):
    mm_count+=1
  if (aa4_donor != aa3_recip):
    mm_count+=1
  return mm_count

# Count number of mismatches between alleles at a given position, considering
# DQA1 and DQB1 combinations
def count_AA_Mismatches_Allele_DQ(self,
                                  allele1_donor,allele2_donor,allele3_donor,allele4_donor,
                                  allele1_recip,allele2_recip,allele3_recip,allele4_recip,
                                  position):
  donor_homoz = 0
  if (allele1_donor == allele2_donor):
    donor_homoz+=1
  if (allele3_donor == allele4_donor):
    donor_homoz+=1
  print ("Number of homozygous donor loci: " + str(donor_homoz))

  aa1_donor = self.getAAposition(allele1_donor,position)
  aa2_donor = self.getAAposition(allele2_donor,position)
  aa3_donor = self.getAAposition(allele3_donor,position)
  aa4_donor = self.getAAposition(allele4_donor,position)
  aa1_recip = self.getAAposition(allele1_recip,position)
  aa2_recip = self.getAAposition(allele2_recip,position)
  aa3_recip = self.getAAposition(allele3_recip,position)
  aa4_recip = self.getAAposition(allele4_recip,position)

  print(aa1_donor)
  print(aa2_donor)
  print(aa3_donor)
  print(aa4_donor)
  print(aa1_recip)
  print(aa2_recip)
  print(aa3_recip)
  print(aa4_recip)
  
  mm_count = self.count_AA_Mismatches_DQ(aa1_donor,aa2_donor,aa3_donor,aa4_donor,aa1_recip,aa2_recip,aa3_recip,aa4_recip)

  if (donor_homoz == 1):
    mm_count-=2
  if (donor_homoz == 2):
    mm_count-=4
  print ("Number of AAMM at position " + str(position) + " : " + str(mm_count))

  if mm_count == 0:
    pass
  elif mm_count >= 1:
    if position in DQ_eplet_positions_all:
      return position
    else:
      pass

# Test DQ AAMM computation, considering DQA1B1 combinations at all positions

allele1_donor = "DQA1*01:01"
allele2_donor = "DQA1*03:01"
allele3_donor = "DQB1*05:01"
allele4_donor = "DQB1*03:02"
allele1_recip = "DQA1*01:02"
allele2_recip = "DQA1*05:01"
allele3_recip = "DQB1*06:04"
allele4_recip = "DQB1*02:01"

position = 13

aa_mm = AAMatch(dbversion=3420)
aam = aa_mm.count_AA_Mismatches_Allele_DQ(allele1_donor,allele2_donor,allele3_donor,allele4_donor,
                                          allele1_recip,allele2_recip,allele3_recip,allele4_recip,position)

print(aam)