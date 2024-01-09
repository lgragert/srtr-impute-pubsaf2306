# Add allele MM to the antigen MM file
import pandas as pd
import numpy as np
import re

# Add the SRTR computed HLA loci A, B, and DR from the TX_KI_decoded.txt
srtr = pd.read_csv("TX_KI_decoded.txt", sep='\t')

SRTR = srtr[['PX_ID', 'DON_A1', 'DON_A2', 'REC_A1', 'REC_A2', 'REC_A_MM_EQUIV_CUR', 'DON_B1', 'DON_B2', 'REC_B1',
             'REC_B2', 'REC_B_MM_EQUIV_CUR', 'DON_DR1', 'DON_DR2', 'REC_DR1', 'REC_DR2',  'REC_DR_MM_EQUIV_CUR']]


# Fix HLA typings into A*01:01 format
def take_spaces_out(SRTR, letter):
    SRTR_new = SRTR
    if letter == "DR":
        SRTR_new[letter] = letter + 'B1*'
    else:
        SRTR_new[letter] = letter + '*'
    column_names = ['DON_' + letter + '1', 'DON_' + letter + '2', 'REC_' + letter + '1', "REC_" + letter + '2']
    for col in column_names:
        SRTR_new[col] = SRTR[col].str.replace('          ', '')
        SRTR_new[col] = SRTR_new[col].str.replace(': ', ':')
        SRTR_new[col] = SRTR_new[col].str.replace(' ', '')

        SRTR_new[col] = SRTR_new[letter] + SRTR_new[col]

    return SRTR_new


SRTR_A = take_spaces_out(SRTR, "A")
SRTR_B = take_spaces_out(SRTR_A, 'B')
SRTR_DR = take_spaces_out(SRTR_B, 'DR')

SRTR_DR = SRTR_DR.drop(['A', 'B', 'DR'], axis=1)


# Compute allele MM by comparing strings (should get 0, 1, 2), column name should be like REC_A_ALLELE_MM
def allele_counting(AG_MM, letter):
    donor = 'DON_' + letter
    recp = 'REC_' + letter

    if letter == 'DQA' or letter == 'DPA' or letter == 'DPB':
        letter = letter + '1'

    AG_MM['REC1_' + letter + '_ALLELE_MM'] = AG_MM[donor + '1'] != AG_MM[recp + '1']
    AG_MM['REC2_' + letter + '_ALLELE_MM'] = AG_MM[donor + '2'] != AG_MM[recp + '2']
    AG_MM['REC_' + letter + '_ALLELE_MM'] = AG_MM[['REC1_' + letter + '_ALLELE_MM', 'REC2_' + letter + '_ALLELE_MM']].sum(axis=1)

    AG_MM = AG_MM.drop(['REC1_' + letter + '_ALLELE_MM', 'REC2_' + letter + '_ALLELE_MM'], axis=1)

    return AG_MM


letters = ['A', 'B', 'C', 'DR', 'DQ', 'DQA', 'DPA', 'DPB']
for num in range(1,11):
    filename = "srtr_antigen_mm_" + str(num) + ".csv"
    ag_file = pd.read_csv(filename)

    # Merge the 3 loci from before with the other loci
    ag_mm = pd.merge(ag_file, SRTR_DR, how='inner', on='PX_ID')

    ag_mm_allele = ag_mm

    for allele in letters:
        ag_mm_allele = allele_counting(ag_mm_allele, allele)

    # Fix the format so that it looks better
    ag_mm_allele = ag_mm_allele[['PX_ID', 'DON_A1', 'DON_A2', 'REC_A1', 'REC_A2', 'REC_A_MM_EQUIV_CUR', 'REC_A_ALLELE_MM',
                                 'DON_B1', 'DON_B2', 'REC_B1', 'REC_B2', 'REC_B_MM_EQUIV_CUR', 'REC_B_ALLELE_MM',
                                 'DON_C1', 'DON_C2', 'REC_C1', 'REC_C2', 'REC_C_MM_EQUIV_CUR', 'REC_C_ALLELE_MM',
                                 'DON_DR1', 'DON_DR2', 'REC_DR1', 'REC_DR2', 'REC_DR_MM_EQUIV_CUR', 'REC_DR_ALLELE_MM',
                                 'DON_DQ1', 'DON_DQ2', 'REC_DQ1', 'REC_DQ2', 'REC_DQ_MM_EQUIV_CUR', 'REC_DQ_ALLELE_MM',
                                 'DON_DQA1', 'DON_DQA2', 'REC_DQA1', 'REC_DQA2', 'REC_DQA1_MM_EQUIV_CUR', 'REC_DQA1_ALLELE_MM',
                                 'DON_DPA1', 'DON_DPA2', 'REC_DPA1', 'REC_DPA2', 'REC_DPA1_MM_EQUIV_CUR', 'REC_DPA1_ALLELE_MM',
                                 'DON_DPB1', 'DON_DPB2', 'REC_DPB1', 'REC_DPB2', 'REC_DPB1_MM_EQUIV_CUR', 'REC_DPB1_ALLELE_MM']]

    ag_mm_allele.to_csv('srtr_ag_allele_mm_' + str(num) + '.csv', header=True, index=False)




