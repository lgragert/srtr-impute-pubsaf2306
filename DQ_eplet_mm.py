import pandas as pd
import re
from aa_matching_msf_genie_MOD import *

# All DQ eplets listed in Eplet Registry
DQ_epReg = pd.read_csv("Eplet_Registry_DQ.txt", sep='\t') 

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

# Make eplet polymorphism computable, listing out all positions in eplet
def polymorph(EpName):
    # Find the first occurrence of a number in the string
    positionInt = re.search(r'\d+', EpName)
    if not positionInt:
        return EpName  # If no number is found, return the original string
    # Extract the starting number and its position (index)
    start_number = int(positionInt.group())
    start_index = positionInt.start()
    end_index = positionInt.end()
    # Initialize the result with the starting number and skip preceding letters
    polym_str = []
    i = start_index
    while i < len(EpName):
        if EpName[i].isdigit():
            # Add the number to the result and skip continuous digits
            if i == start_index:
                polym_str.append(str(start_number))
                start_number += 1
                i = end_index - 1
        else:
            # Add the letter to the result
            if EpName[i].isalpha() and EpName[i-1].isdigit():
              polym_str.append(EpName[i])
            elif EpName[i].isalpha() and EpName[i-1].isalpha():
              polym_str.append(str(start_number))
              start_number += 1
              polym_str.append(EpName[i])
            elif not EpName[i].isalnum():
              polym_str.append(EpName[i])
            elif EpName[i].isalpha() and not EpName[i-1].isalnum():
              polym_str.append(EpName[i])
        i += 1
    return ''.join(polym_str)
# Apply the function to the 'Eplet_Name' column to create a new column
DQA1_eplets['Polymorphism'] = DQA1_eplets['Eplet_Name'].apply(polymorph)
DQB1_eplets['Polymorphism'] = DQB1_eplets['Eplet_Name'].apply(polymorph)

# Create dictionary for all DQ eplets in eplet registry
DQ_eplet_positions_all = {}  # all positions in eplet registry
DQ_eplet_positions_calculator = {}  # all positions used by eplet MM calculator
for eplet_polymorphism in eplet_DQ['Polymorphic']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_all[position] = 1
    # Check if is in hlaR reference data
    # if statement then DQ_eplet_positions_calculator[position] = 1
# Created a dictionary of eplet positions listed in Eplet_Registry_DQ.txt
# Key of AA position
# Value of 1 if it has been seen
print(DQ_eplet_positions_all)

# Create dictionary for DQ eplets used by calculator
DQ_eplet_positions_calculator_A = {} # A positions used by eplet MM calculator
DQ_eplet_positions_calculator_B = {} # B positions used by eplet MM calculator

for eplet_polymorphism in DQA1_eplets['Polymorphism']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_calculator_A[position] = 1

for eplet_polymorphism in DQB1_eplets['Polymorphism']:
  eplet_positions = re.findall(r'\d+', eplet_polymorphism)
  for position in eplet_positions:
    DQ_eplet_positions_calculator_B[position] = 1

# Combine both dictionaries
DQ_eplet_positions_calculator_str = DQ_eplet_positions_calculator_A | DQ_eplet_positions_calculator_B
# Convert dictionary key strings into integers
DQ_eplet_positions_calculator = {int(key): value for key, value 
                                 in DQ_eplet_positions_calculator_str.items()}

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

position = 13 # 0 AAMMs, returns None
aa_mm = AAMatch(dbversion=3420)
pos13 = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                          DQB1_1_donor,DQB1_2_donor,
                                          DQA1_1_recip,DQA1_2_recip,
                                          DQB1_1_recip,DQB1_2_recip,position)
print(pos13)

position = 75 # 4 AAMMs, is eplet, returns None
pos75 = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                          DQB1_1_donor,DQB1_2_donor,
                                          DQA1_1_recip,DQA1_2_recip,
                                          DQB1_1_recip,DQB1_2_recip,position)
print(pos75)

position = 11 # 2 AAMMs, non-eplet, returns position 11
pos11 = aa_mm.count_AA_Mismatches_Allele_DQ(DQA1_1_donor,DQA1_2_donor,
                                          DQB1_1_donor,DQB1_2_donor,
                                          DQA1_1_recip,DQA1_2_recip,
                                          DQB1_1_recip,DQB1_2_recip,position)
print(pos11)

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
