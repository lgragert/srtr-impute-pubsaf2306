import pandas as pd
import re
from aa_matching_msf_genie_MOD import *

# Read in DQA1 and DQB1 eplets from Eplet Registry (Original)
DQ_epReg1_full = pd.read_csv("Eplet_Registry_DQ.txt", sep='\t')
DQ_epReg1 = DQ_epReg1_full[['Eplet', 'Polymorphic']]
# Read in DQA1 and DQB1 eplets from Eplet Registry (Updated)
DQ_epReg2_full = pd.read_csv("Eplet_Registry_DQ_Updated.csv")
DQ_epReg2 = DQ_epReg2_full[['Name', 'Description']]
DQ_epReg2.columns = DQ_epReg1.columns
# Read in Interlocus eplets from Eplet Registry
Inter_epReg_full = pd.read_csv("Eplet_Registry_Interlocus.txt", sep='\t',
                               encoding_errors='ignore')
Inter_epReg = Inter_epReg_full[['Name', 'Description']]
Inter_epReg.columns = DQ_epReg1.columns
# Combine all DQ eplets from Eplet Registry
AllDQ_epReg = pd.concat([DQ_epReg1, DQ_epReg2, Inter_epReg], ignore_index=True)

# Read in DQA1 eplets from LarsenLab hlaR
DQ_hlaRA = pd.read_csv("MHC_II_eplet_A_v3.csv")
# Filter and remove non-DQ columns
DPA1_columns = DQ_hlaRA.filter(regex='^DP').columns
DQ_hlaRA = DQ_hlaRA.drop(columns=DPA1_columns)
# Repeat for DQB1 eplets
DQ_hlaRB = pd.read_csv("MHC_II_eplet_B_v3.csv")
DPB1_DRB1345_columns = DQ_hlaRB.filter(regex='^(DP|DR)').columns
DQ_hlaRB = DQ_hlaRB.drop(columns=DPB1_DRB1345_columns)

# Identify DQA1 and DQB1 eplets (positions)
# Initialize an empty list to store non-NaN values
DQA1_eplet_names = []
for col in DQ_hlaRA.columns[2:]:
  # Extract non-NaN values from the column and add them to the list
  DQA1_eplet_names.extend(DQ_hlaRA[col].dropna().unique())
# To sort positions in order (numerically)
def custom_sort(value):
  # Extract numeric part of the value
  numeric_part = int(''.join(filter(str.isdigit, value)))
  return numeric_part
# Remove repeated values and sort them in numerical order
# using the custom sorting function
DQA1_eplet_names_sorted = sorted(set(DQA1_eplet_names), key=custom_sort)
# Create a new DataFrame with a single column containing the sorted values
DQA1_eplets = pd.DataFrame(DQA1_eplet_names_sorted, columns=['Eplet_Name'])
# Repeat for DQB1 eplets
DQB1_eplet_names = []
for col in DQ_hlaRB.columns[2:]:
  DQB1_eplet_names.extend(DQ_hlaRB[col].dropna().unique())
DQB1_eplet_names_sorted = sorted(set(DQB1_eplet_names), key=custom_sort)
DQB1_eplets = pd.DataFrame(DQB1_eplet_names_sorted, columns=['Eplet_Name'])
DQB1_eplets.at[32, 'Eplet_Name'] = 'rq70RK/R'
# Combine DQA1 and DQB1 eplet names
DQA1B1_hlaR = pd.concat([DQA1_eplets, DQB1_eplets], ignore_index=True)

# Make eplet polymorphism computable, listing out all positions in eplet
# Convert Eplet_Name column to uppercase for case-insensitive comparison
DQA1B1_hlaR['Eplet_Name'] = DQA1B1_hlaR['Eplet_Name'].str.upper()

DQA1B1_hlaR_merged = pd.merge(DQA1B1_hlaR, AllDQ_epReg, how='left',
                              left_on='Eplet_Name', right_on='Eplet')

# hlaR eplets that are not listed in Eplet Registry
# Input single AA polymorphic strings
hlaR_polymorphic = [('PQ34Q', '34Q'), ('69L', '69L'), ('69T', '69T'),
                    ('129Q', '129Q'), ('130S', '130S'), ('66E', '66E'),
                    ('67V', '67V'), ('67I', '67I'), ('67D', '67D'),
                    ('70R', '70R'), ('70E', '70E'), ('70G', '70G'),
                    ('71D', '71D'), ('71T', '71T'), ('71A', '71A'),
                    ('71K', '71K'), ('86G', '86G')]
for eplet_identifier, pmrph_input in hlaR_polymorphic:
  DQA1B1_hlaR_merged.loc[DQA1B1_hlaR_merged['Eplet_Name'] == eplet_identifier,
                         'Polymorphic'] = pmrph_input

DQA1B1_hlaR_merged_filt = DQA1B1_hlaR_merged.dropna(subset=['Polymorphic'])

# Create dictionary for all DQ eplets in eplet registry
DQ_eplet_positions_all_str = {}  # all positions in eplet registry
DQ_eplet_positions_calculator = {}  # all positions used by eplet MM calculator
for eplet_polymorphism in AllDQ_epReg['Polymorphic']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_all_str[position] = 1
    # Check if is in hlaR reference data
    # if statement then DQ_eplet_positions_calculator[position] = 1
# Created a dictionary of eplet positions listed in Eplet_Registry_DQ.txt
# Key of AA position
# Value of 1 if it has been seen
print(DQ_eplet_positions_all)

# Create dictionary for DQ eplets used by calculator
DQ_eplet_positions_calculator_str = {} # all positions used by eplet MM calculator

for eplet_polymorphism in DQA1B1_hlaR_merged_filt['Polymorphic']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_calculator_str[position] = 1

# Find keys in Dict(DQ_eplet_positions_calculator_str) that are
# NOT present in Dict(DQ_eplet_positions_all_str)
hlaR_only_eps_str = [key for key in DQ_eplet_positions_calculator_str
                     if key not in DQ_eplet_positions_all_str]
print(hlaR_only_eps_str)

# Convert dictionary key strings into integers
DQ_eplet_positions_calculator = {int(key): value for key, value 
                                 in DQ_eplet_positions_calculator_str.items()}
hlaR_only_eps = [int(item) for item in hlaR_only_eps_str]

# Test DQ AAMM computation, considering DQA1B1 combinations at all positions
# Rename these to DQA1_1_donor to agree with input variables when you call the function below
# input variables do not need to agree with the variable names in the function definition
DQA1_1_donor = "DQA1*01:01"
DQA1_2_donor = "DQA1*03:01"
DQB1_1_donor = "DQB1*05:01"
DQB1_2_donor = "DQB1*03:02"
DQA1_1_recip = "DQA1*01:02"
DQA1_2_recip = "DQA1*05:01"
DQB1_1_recip = "DQB1*06:04"
DQB1_2_recip = "DQB1*02:01"

position = 14 # 1 AAMM, returns position 14
aa_mm = AAMatch(dbversion=3420)
pos14 = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                            DQB1_1_donor,DQB1_2_donor,
                                            DQA1_1_recip,DQA1_2_recip,
                                            DQB1_1_recip,DQB1_2_recip,position)
print(pos14)

position = 34 # is eplet or non-eplet, returns position 34
pos34 = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                            DQB1_1_donor,DQB1_2_donor,
                                            DQA1_1_recip,DQA1_2_recip,
                                            DQB1_1_recip,DQB1_2_recip,position)
print(pos34)

# Loop through all eplet positions
# Initialize non-eplet MM positions
nonEp = []
for position in range(1, 201):
  allpos = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                               DQB1_1_donor,DQB1_2_donor,
                                               DQA1_1_recip,DQA1_2_recip,
                                               DQB1_1_recip,DQB1_2_recip,position)
  nonEp.append((position, allpos))

nonEp_df = pd.DataFrame(nonEp, columns=['Position', 'NonEplet_Position'])
nonEp_df_filt = nonEp_df.dropna(subset=['NonEplet_Position'])
print(nonEp_df_filt)
