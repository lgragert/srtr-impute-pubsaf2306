import pandas as pd
import re
from aa_matching_msf_genie_MOD import *

eplet_DQ = pd.read_csv("Eplet_Registry_DQ.txt", sep='\t') 

DQ_eplet_positions_all = {}  # all positions in eplet registry
DQ_eplet_positions_calculator = {}  # all positions used by eplet MM calculator

for eplet_polymorphism in eplet_DQ['Polymorphic']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_all[position] = 1
    # TODO - check if is in hlaR reference data
    # if statement then DQ_eplet_positions_calculator[position] = 1

# Created a dictionary of eplet positions listed in Eplet_Registry_DQ.txt
# Key of AA position
# Value of 1 if it has been seen
print(DQ_eplet_positions_all)

# Test DQ AAMM computation, considering DQA1B1 combinations at all positions

# TODO - rename these to DQA1_1_donor to agree with input variables when you call the function below
# input variables do not need to agree with the variable names in the function definition
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
aam = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                          DQB1_1_donor,DQB1_2_donor,
                                          DQA1_1_recip,DQA1_2_recip,
                                          DQB1_1_recip,DQB1_2_recip,position)

print(aam)

position = 75
aam2 = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                          DQB1_1_donor,DQB1_2_donor,
                                          DQA1_1_recip,DQA1_2_recip,
                                          DQB1_1_recip,DQB1_2_recip,position)
print(aam2)

#NameError: name 'DQ_eplet_positions_all' is not defined

for i in range(1, 201):
  count_AA_Mismatches_Allele_DQ(self,
                                  DQA1_1_donor,DQA1_2_donor,
                                  DQB1_1_donor,DQB1_2_donor,
                                  DQA1_1_recip,DQA1_2_recip,
                                  DQB1_1_recip,DQB1_2_recip,
                                  i)
