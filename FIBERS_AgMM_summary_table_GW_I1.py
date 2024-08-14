import pandas as pd
import os
import sys

# Script to generate summary stats for antigen MM and amino acid MM

risk_group = sys.argv[1] # OVERALL, FIBERS_LOW, FIBERS_HIGH

# SRTR imputation replicate 1
SRTR_imputation_replicate1_filename = "./SRTR_AA_MM_matrix_grffail_replicates_2022-04-03/SRTR_AA_MM_matrix_grffail_1.txt" 
SRTR = pd.read_csv(SRTR_imputation_replicate1_filename,sep='\t')

# Alternative imputation realization 9
# SRTR_imputation_replicate9_filename = "SRTR_AA_MM_matrix_grffail_9.txt" 
# SRTR_R9 = pd.read_csv(SRTR_imputation_replicate9_filename,sep='\t')


# subset antigen MM and amino acid MM columns from data frame

# FIBERS Bin Imputation 9 AAs - What we are calling the "Top FIBERS bin"
# "MM_DRB1_13","MM_DRB1_26","MM_DQB1_30","MM_DQB1_55"

# FIBERS Bin Imputation 1 AAs - A different larger set of AAs
# "MM_DRB1_11","MM_DRB1_26","MM_DRB1_37","MM_DQB1_53","MM_DRB1_77","MM_DQB1_89"

SRTR = SRTR[["PX_ID", "REC_DR_MM_EQUIV_CUR", "REC_DQ_MM_EQUIV_CUR",
            "MM_DRB1_13","MM_DRB1_26","MM_DQB1_30","MM_DQB1_55",
            "MM_DRB1_11","MM_DRB1_37","MM_DQB1_53","MM_DRB1_77","MM_DQB1_89"
            ]]


# overall count
kidney_tx_pair_count = len(SRTR)
print ("SRTR Mismatch Summary Statistics:\n")
print ("Overall Count: " + str(kidney_tx_pair_count) + "\n")

# FIBERS Risk Categories

print ("SRTR Deceased Kidney Donors")
print ("Imputation Replicate 1")
print ("FIBERS - AAMM Positions in Bin: MM_DRB1_11,MM_DRB1_26,MM_DRB1_37,MM_DQB1_53,MM_DRB1_77,MM_DQB1_89")

SRTR_Low_Risk = SRTR[(SRTR['MM_DRB1_11'] == 0) &
                     (SRTR['MM_DRB1_26'] == 0) &
                     (SRTR['MM_DRB1_37'] == 0) &
                     (SRTR['MM_DRB1_77'] == 0) &
                     (SRTR['MM_DQB1_53'] == 0) &
                     (SRTR['MM_DQB1_89'] == 0)]

FIBERS_Low_Risk_count = len(SRTR_Low_Risk)

print ("FIBERS Low Risk Count: " + str(FIBERS_Low_Risk_count))


SRTR_High_Risk = SRTR[(SRTR['MM_DRB1_11'] >= 1) |
                      (SRTR['MM_DRB1_26'] >= 1) |
                      (SRTR['MM_DRB1_37'] >= 1) |
                      (SRTR['MM_DRB1_77'] >= 1) |
                      (SRTR['MM_DQB1_53'] >= 1) |
                      (SRTR['MM_DQB1_89'] >= 1) ]

FIBERS_High_Risk_count = len(SRTR_High_Risk)

print ("FIBERS High Risk Count: " + str(FIBERS_High_Risk_count))

# select appropriate data frame for risk group
if (risk_group == "OVERALL"):
    print ("\nSTATS FOR OVERALL DECEASED DONOR KIDNEY TX PAIRS\n")
elif (risk_group == "FIBERS_LOW"):
    print ("\nSTATS WITHIN FIBERS LOW RISK GROUP (0 AAMM in bin)\n")
    SRTR = SRTR_Low_Risk
    kidney_tx_pair_count = len(SRTR)
elif (risk_group == "FIBERS_HIGH"):
    print ("\nSTATS WITHIN FIBERS HIGH RISK GROUP (1+ AAMM in bin)\n")
    SRTR = SRTR_High_Risk
    kidney_tx_pair_count = len(SRTR)
else:
    print ("Error - Invalid Risk Group: " + risk_group)
    exit()


# Convert AgMM to integers and handle NAs due to missing DQ typing
# 99 - Not tested per SRTR data dictionary
SRTR['REC_DR_MM_EQUIV_CUR'] = SRTR['REC_DR_MM_EQUIV_CUR'].fillna(99).astype(int)
SRTR['REC_DQ_MM_EQUIV_CUR'] = SRTR['REC_DQ_MM_EQUIV_CUR'].fillna(99).astype(int)
SRTR = SRTR.astype({'REC_DR_MM_EQUIV_CUR':'int','REC_DQ_MM_EQUIV_CUR':'int'})





# Overall Antigen Mismatch Counts

print ("\nAntigen Mismatch Counts\n")

antigen_DR_NA_count = len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 99])
antigen_DQ_NA_count = len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 99])
antigen_DR_NA_fract = antigen_DR_NA_count / kidney_tx_pair_count
antigen_DQ_NA_fract = antigen_DQ_NA_count / kidney_tx_pair_count

print ("Missing-DR AgMM Count: " + str(antigen_DR_NA_count) + " (" + 
            format(antigen_DR_NA_fract*100, '.2f') + "%)")
print ("Missing-DQ AgMM Count: " + str(antigen_DQ_NA_count) + " (" + 
            format(antigen_DQ_NA_fract*100, '.2f') + "%)")

antigen_DR_0MM_count = len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 0])
antigen_DQ_0MM_count = len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 0])
antigen_DR_0MM_fract = antigen_DR_0MM_count / kidney_tx_pair_count
antigen_DQ_0MM_fract = antigen_DQ_0MM_count / kidney_tx_pair_count

print ("0-DR AgMM Count: " + str(antigen_DR_0MM_count) + " (" + 
            format(antigen_DR_0MM_fract*100, '.2f') + "%)")
print ("0-DQ AgMM Count: " + str(antigen_DQ_0MM_count) + " (" + 
            format(antigen_DQ_0MM_fract*100, '.2f') + "%)")


antigen_DR_1MM_count = len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 1])
antigen_DQ_1MM_count = len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 1])
antigen_DR_1MM_fract = antigen_DR_1MM_count / kidney_tx_pair_count
antigen_DQ_1MM_fract = antigen_DQ_1MM_count / kidney_tx_pair_count

print ("1-DR AgMM Count: " + str(antigen_DR_1MM_count) + " (" + 
            format(antigen_DR_1MM_fract*100, '.2f') + "%)")
print ("1-DQ AgMM Count: " + str(antigen_DQ_1MM_count) + " (" + 
            format(antigen_DQ_1MM_fract*100, '.2f') + "%)")

antigen_DR_2MM_count = len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 2])
antigen_DQ_2MM_count = len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 2])
antigen_DR_2MM_fract = antigen_DR_2MM_count / kidney_tx_pair_count
antigen_DQ_2MM_fract = antigen_DQ_2MM_count / kidney_tx_pair_count

print ("2-DR AgMM Count: " + str(antigen_DR_2MM_count) + " (" + 
            format(antigen_DR_2MM_fract*100, '.2f') + "%)")
print ("2-DQ AgMM Count: " + str(antigen_DQ_2MM_count) + " (" + 
            format(antigen_DQ_2MM_fract*100, '.2f') + "%)")

antigen_DR_0MM_DQ_0MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 0)])
antigen_DR_0MM_DQ_0MM_fract = antigen_DR_0MM_DQ_0MM_count / kidney_tx_pair_count

print ("0-DR and 0-DQ AgMM Count: " + str(antigen_DR_0MM_DQ_0MM_count) + " (" + 
            format(antigen_DR_0MM_DQ_0MM_fract*100, '.2f') + "%)")

antigen_DR_0MM_DQ_1MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 1)])
antigen_DR_0MM_DQ_1MM_fract = antigen_DR_0MM_DQ_1MM_count / kidney_tx_pair_count

print ("0-DR and 1-DQ AgMM Count: " + str(antigen_DR_0MM_DQ_1MM_count) + " (" + 
            format(antigen_DR_0MM_DQ_1MM_fract*100, '.2f') + "%)")

antigen_DR_1MM_DQ_1MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 1)])
antigen_DR_1MM_DQ_1MM_fract = antigen_DR_1MM_DQ_1MM_count / kidney_tx_pair_count

print ("1-DR and 1-DQ AgMM Count: " + str(antigen_DR_1MM_DQ_1MM_count) + " (" + 
            format(antigen_DR_1MM_DQ_1MM_fract*100, '.2f') + "%)")

antigen_DR_0MM_DQ_2MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 2)])
antigen_DR_0MM_DQ_2MM_fract = antigen_DR_0MM_DQ_2MM_count / kidney_tx_pair_count

print ("0-DR and 2-DQ AgMM Count: " + str(antigen_DR_0MM_DQ_2MM_count) + " (" + 
            format(antigen_DR_0MM_DQ_2MM_fract*100, '.2f') + "%)")

antigen_DR_2MM_DQ_0MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 0)])
antigen_DR_2MM_DQ_0MM_fract = antigen_DR_2MM_DQ_0MM_count / kidney_tx_pair_count

print ("2-DR and 0-DQ AgMM Count: " + str(antigen_DR_2MM_DQ_0MM_count) + " (" + 
            format(antigen_DR_2MM_DQ_0MM_fract*100, '.2f') + "%)")

antigen_DR_1MM_DQ_2MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 2)])
antigen_DR_1MM_DQ_2MM_fract = antigen_DR_1MM_DQ_2MM_count / kidney_tx_pair_count

print ("1-DR and 2-DQ AgMM Count: " + str(antigen_DR_1MM_DQ_2MM_count) + " (" + 
            format(antigen_DR_1MM_DQ_2MM_fract*100, '.2f') + "%)")

antigen_DR_2MM_DQ_1MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 1)])
antigen_DR_2MM_DQ_1MM_fract = antigen_DR_2MM_DQ_1MM_count / kidney_tx_pair_count

print ("2-DR and 1-DQ AgMM Count: " + str(antigen_DR_2MM_DQ_1MM_count) + " (" + 
            format(antigen_DR_2MM_DQ_1MM_fract*100, '.2f') + "%)")

antigen_DR_2MM_DQ_2MM_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 2)])
antigen_DR_2MM_DQ_2MM_fract = antigen_DR_2MM_DQ_2MM_count / kidney_tx_pair_count

print ("2-DR and 2-DQ AgMM Count: " + str(antigen_DR_2MM_DQ_2MM_count) + " (" + 
            format(antigen_DR_2MM_DQ_2MM_fract*100, '.2f') + "%)")


antigen_DR_0MM_DQ_NA_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 99)])
antigen_DR_0MM_DQ_NA_fract = antigen_DR_0MM_DQ_NA_count / kidney_tx_pair_count

print ("0-DR and Missing-DQ AgMM Count: " + str(antigen_DR_0MM_DQ_NA_count) + " (" + 
            format(antigen_DR_0MM_DQ_NA_fract*100, '.2f') + "%)")

antigen_DR_1MM_DQ_NA_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 99)])
antigen_DR_1MM_DQ_NA_fract = antigen_DR_1MM_DQ_NA_count / kidney_tx_pair_count

print ("1-DR and Missing-DQ AgMM Count: " + str(antigen_DR_1MM_DQ_NA_count) + " (" + 
            format(antigen_DR_1MM_DQ_NA_fract*100, '.2f') + "%)")


antigen_DR_2MM_DQ_NA_count = len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 99)])
antigen_DR_2MM_DQ_NA_fract = antigen_DR_2MM_DQ_NA_count / kidney_tx_pair_count

print ("2-DR and Missing-DQ AgMM Count: " + str(antigen_DR_2MM_DQ_NA_count) + " (" + 
            format(antigen_DR_2MM_DQ_NA_fract*100, '.2f') + "%)")


# VERBOSE TABLES

print ("\nVERBOSE TABLES FROM EXPLORATORY DATA ANALYSIS\n")

# Antigen MM summary stats

print ("\nAntigen Mismatch Summary value_counts\n")

print (SRTR['REC_DR_MM_EQUIV_CUR'].value_counts(dropna=False).sort_index())
print (SRTR['REC_DQ_MM_EQUIV_CUR'].value_counts(dropna=False).sort_index())

# Amino Acid MM summary stats
print ("\nAmino Acid MM Summary value_counts\n")

# FIBERS Bin 9 Amino Acids
print (SRTR['MM_DRB1_13'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DRB1_26'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DQB1_30'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DQB1_55'].value_counts(dropna=False).sort_index())


# FIBERS Bin 1
print ("\nAmino Acid MM Summary value_counts - FIBERS Top Bin Replicate 1 \n")

print (SRTR['MM_DRB1_11'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DRB1_26'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DRB1_37'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DRB1_77'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DQB1_53'].value_counts(dropna=False).sort_index())
print (SRTR['MM_DQB1_89'].value_counts(dropna=False).sort_index())
