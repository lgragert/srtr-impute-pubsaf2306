#!/usr/bin/env python
# Converting SRTR data to PRICHE format

import pandas as pd

# Read in the first replicate
SRTR = pd.read_csv('./SRTR_AA_MM_9loc_matrix_1.txt.gz', sep='\t')

SRTR_pirche_info = SRTR[['PX_ID', 'HAPPAIR_RECIP', 'HAPPAIR_DONOR', 'CAN_RACE', 'DON_RACE']]

# The PIRCHE format is recip, happair endl, donor, happair, endl, endl

# Create a donor DataFrame
donor = SRTR_pirche_info[['PX_ID', 'HAPPAIR_DONOR', 'DON_RACE']].copy()
donor['PX_ID'] = 'D' + donor['PX_ID'].astype(str)
donor = donor.rename(columns={'PX_ID': 'ID', 'HAPPAIR_DONOR': 'HAPPAIR', 'DON_RACE': 'RACE'})

# Create recipient DataFrame
recip = SRTR_pirche_info[['PX_ID', 'HAPPAIR_RECIP', 'CAN_RACE']].copy()
recip['PX_ID'] = 'R' + recip['PX_ID'].astype(str)
recip = recip.rename(columns={'PX_ID': 'ID', 'HAPPAIR_RECIP': 'HAPPAIR', 'CAN_RACE': 'RACE'})


# Concatenate donor, recipient, and blank lines efficiently by sorting the index
recip_donor_df = pd.concat([recip, donor]).sort_index(kind='stable')

# Create a blank line DataFrame that is the same size as the recip
blank_rows = pd.DataFrame({'ID': [''] * len(recip), 'HAPPAIR': [''] * len(recip), 'RACE': [''] * len(recip)})

pirche_dataframe = pd.concat([recip_donor_df, blank_rows], axis=0).sort_index(kind='stable').reset_index(drop=True)


# Expand happair, so there are commas and not glstring format
pirche_dataframe[['HAPPAIR_1', 'HAPPAIR_2']] = pirche_dataframe.HAPPAIR.str.split('+', expand=True)
pirche_dataframe[["A_1", "C_1", "B_1", "DRB345_1", "DRB1_1", "DQA1_1", "DQB1_1", "DPA1_1", "DPB1_1"]] = (
    pirche_dataframe.HAPPAIR_1.str.split('~', expand=True))
pirche_dataframe[["A_2", "C_2", "B_2", "DRB345_2", "DRB1_2", "DQA1_2", "DQB1_2", "DPA1_2", "DPB1_2"]] = (
    pirche_dataframe.HAPPAIR_2.str.split('~', expand=True))

# drop HAPPAIR columns for 9-locus version
PIRCHE_9LOC = pirche_dataframe
PIRCHE_9LOC = PIRCHE_9LOC.drop(columns=['HAPPAIR', 'HAPPAIR_1', 'HAPPAIR_2'])
PIRCHE_9LOC = PIRCHE_9LOC[:-1]  # need to drop the last set of commas
print(PIRCHE_9LOC.head())

# Go from DataFrame to PIRCHE format CSV
PIRCHE_9LOC.to_csv('PIRCHE_SRTR_9loc_1.csv', header=False, index=False)
