import pandas as pd
from aa_matching_msf_genie_MOD import *

# test DR amino acid mismatch computation, considering both DRB1 and DRB3/4/5 loci as homologous

allele1_donor = "DRB1*07:01"
allele2_donor = "DRB1*15:01"
allele3_donor = "DRB4*01:03"
allele4_donor = "DRB5*01:01"
allele1_recip = "DRB1*11:01"
allele2_recip = "DRB1*15:03"
allele3_recip = "DRB3*02:02"
allele4_recip = "DRB4*01:01"

position = 13 # most variable position in DRB1

aa_mm = AAMatch(dbversion=3420)
aam = aa_mm.count_AA_Mismatches_Allele_DR(allele1_donor,allele2_donor,allele3_donor,allele4_donor,allele1_recip,allele2_recip,allele3_recip,allele4_recip,position)

print(aam)

# DRB1345 test cases

testinput = pd.read_csv('SRTR_AA_MM_9loc_matrix_1_EXTR.txt', delimiter='\t')

testinput['DRB1345_P13'] = testinput.apply(lambda row: aa_mm.count_AA_Mismatches_Allele_DR(row['DONOR_DRB1_1'], row['DONOR_DRB1_2'],
                                                                                           row['DONOR_DRB345_1'], row['DONOR_DRB345_2'],
                                                                                           row['RECIP_DRB1_1'], row['RECIP_DRB1_2'],
                                                                                           row['RECIP_DRB345_1'], row['RECIP_DRB345_2'],
                                                                                           position), axis=1)

testinput

def summary_statistics(testinput):
  # Convert tuples to strings in [DRB1345_P13] column
  # Extract number of mismatches from [DRB1345_P13] column
  testinput['DRB1345_P13_COUNT'] = testinput['DRB1345_P13'].apply(lambda x: int(str(x).split(',')[1].strip(')')))

  # Summary statistics (overall)
  SS_overall = testinput['DRB1345_P13_COUNT'].value_counts().sort_index()

  # Summary statistics (by recipient population category in [CAN_RACE] column)
  SS_rpop = testinput.groupby('CAN_RACE')['DRB1345_P13_COUNT'].value_counts().unstack(fill_value=0).sort_index()

  # Print overall summary
  print("Overall Summary Statistics:")
  print(SS_overall)

  # Print summary by recipient population category
  print("\nSummary Statistics by Recipient Population Category:")
  print(SS_rpop)

# Call the function with the DataFrame
summary_statistics(testinput)

# testinput.to_csv('SRTR_AA_MM_9loc_matrix_1_EXTR_DRBoutput.txt', index=False, sep='\t')
