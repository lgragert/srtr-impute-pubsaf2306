import pandas as pd
import os
import sys
from openpyxl.styles import Alignment, Font
from openpyxl import load_workbook
import json

# Script to generate summary stats for antigen MM and amino acid MM


# Format for mismatch counts, this is what will appear on the excel spreadsheet
def count_format(count, total_sum):
    # Puts count into a str format plus with a percentage out of the total number of counts
    ag_count_perc = str(count) + " (" + str(format(100 * count / total_sum, ".2f")) + "%)"
    return ag_count_perc


# Overall Antigen Mismatch Counts Function
def antigen_count(SRTR, replicate, total_count):

    antigen_DR_NA_perc = count_format(len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 99]), total_count)  # Missing DR count
    antigen_DQ_NA_perc = count_format(len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 99]), total_count)  # Missing DQ count

    antigen_ABDR_0MM_perc = count_format(len(SRTR[(SRTR['REC_A_MM_EQUIV_CUR'] == 0) &
                                                  (SRTR['REC_B_MM_EQUIV_CUR'] == 0) &
                                                  (SRTR['REC_DR_MM_EQUIV_CUR'] == 0)]), total_count)  # 0MM at -A,-B,-DR

    antigen_DR_0MM_perc = count_format(len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 0]), total_count)  # 0MM DR
    antigen_DQ_0MM_perc = count_format(len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 0]), total_count)  # 0MM DQ

    antigen_DR_1MM_perc = count_format(len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 1]), total_count)  # 1MM DR
    antigen_DQ_1MM_perc = count_format(len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 1]), total_count)  # 1MM DQ

    antigen_DR_2MM_perc = count_format(len(SRTR[SRTR['REC_DR_MM_EQUIV_CUR'] == 2]), total_count)  # 2MM DR
    antigen_DQ_2MM_perc = count_format(len(SRTR[SRTR['REC_DQ_MM_EQUIV_CUR'] == 2]), total_count)  # 2MM DQ

    antigen_DR_0MM_DQ_0MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 0)]), total_count)  # 0MM DR, DQ

    antigen_DR_0MM_DQ_1MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 1)]), total_count)  # 0MM DR, 1MM DQ

    antigen_DR_1MM_DQ_0MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 0)]), total_count)  # 1MM DR, 0MM DQ

    antigen_DR_1MM_DQ_1MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 1)]), total_count)  # 1MM DR, DQ

    antigen_DR_0MM_DQ_2MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 2)]), total_count)  # 0MM DR, 2MM DQ

    antigen_DR_2MM_DQ_0MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 0)]), total_count)  # 2MM DR, 0MM DQ

    antigen_DR_1MM_DQ_2MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 2)]), total_count)  # 1MM DR, 2MM DQ

    antigen_DR_2MM_DQ_1MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 1)]), total_count)  # 2MM DR, 1MM DQ

    antigen_DR_2MM_DQ_2MM_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                       (SRTR['REC_DQ_MM_EQUIV_CUR'] == 2)]), total_count)  # 2MM DR, DQ

    antigen_DR_0MM_DQ_NA_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                      (SRTR['REC_DQ_MM_EQUIV_CUR'] == 99)]), total_count)  # 0MM DR, Missing DQ

    antigen_DR_1MM_DQ_NA_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                      (SRTR['REC_DQ_MM_EQUIV_CUR'] == 99)]), total_count)  # 1MM DR, Missing DQ

    antigen_DR_2MM_DQ_NA_perc = count_format(len(SRTR[(SRTR['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                      (SRTR['REC_DQ_MM_EQUIV_CUR'] == 99)]), total_count)  # 2MM DR, Missing DQ

    antigen_dict = {"Missing-DR AgMM Count": antigen_DR_NA_perc, "Missing-DQ AgMM Count": antigen_DQ_NA_perc,
                    "0-ABDR AgMM Count": antigen_ABDR_0MM_perc, "0-DR AgMM Count": antigen_DR_0MM_perc,
                    "0-DQ AgMM Count": antigen_DQ_0MM_perc, "1-DR AgMM Count": antigen_DR_1MM_perc,
                    "1-DQ AgMM Count": antigen_DQ_1MM_perc, "2-DR AgMM Count": antigen_DR_2MM_perc,
                    "2-DQ AgMM Count": antigen_DQ_2MM_perc, "0-DR and 0-DQ AgMM Count": antigen_DR_0MM_DQ_0MM_perc,
                    "0-DR and 1-DQ AgMM Count": antigen_DR_0MM_DQ_1MM_perc, "1-DR and 0-DQ AgMM Count": antigen_DR_1MM_DQ_0MM_perc,
                    "1-DR and 1-DQ AgMM Count": antigen_DR_1MM_DQ_1MM_perc, "0-DR and 2-DQ AgMM Count": antigen_DR_0MM_DQ_2MM_perc,
                    "2-DR and 0-DQ AgMM Count": antigen_DR_2MM_DQ_0MM_perc, "1-DR and 2-DQ AgMM Count": antigen_DR_1MM_DQ_2MM_perc,
                    "2-DR and 1-DQ AgMM Count": antigen_DR_2MM_DQ_1MM_perc, "2-DR and 2-DQ AgMM Count": antigen_DR_2MM_DQ_2MM_perc,
                    "0-DR and Missing-DQ AgMM Count": antigen_DR_0MM_DQ_NA_perc, "1-DR and Missing-DQ AgMM Count": antigen_DR_1MM_DQ_NA_perc,
                    "2-DR and Missing-DQ AgMM Count": antigen_DR_2MM_DQ_NA_perc}

    antigen_dataframe = pd.DataFrame(antigen_dict.values, index=list(antigen_dict.keys()), columns=[replicate])

    return antigen_dataframe


# Initialize DataFrames which will be each Excel sheet
headings = ["Replicate_1", "Replicate_2", "Replicate_3", "Replicate_4", "Replicate_5", "Replicate_6", "Replicate_7",
            "Replicate_8", "Replicate_9", "Replicate_10"]
OVERALL = pd.DataFrame(columns=headings)
FIBERS_HIGH = pd.DataFrame(columns=headings)
FIBERS_LOW = pd.DataFrame(columns=headings)

# Get JSON file of all the FIBERS bins with the AAMM that we want
with open('amino_acids.json', 'r') as file:
    config = json.load(file)

FIBERS_AAMM_bin = config['FIBERS_BINS']['Bin_9']  # Only BIN_9 for now
aamm_in_header = ["MM_" + AAMM for AAMM in FIBERS_AAMM_bin]  # Get the column name in matrix files for AAMM
impute_realization = range(1, 11)                 # 1 thru 10 for each realization/replicate

pops = ["OVERALL", "AFA", "CAU", "HIS", "HPI", "NAM", "ASI", "Multi-Racial"]

for pop in pops:
    print("Working on population: ", pop)
    for replicate in impute_realization:
        print("Working on replicate: ", replicate)
        # SRTR imputation replicate
        SRTR_imputation_replicate_filename = "SRTR_AA_MM_9loc_matrix_" + str(replicate) + ".txt.gz"
        SRTR = pd.read_csv(SRTR_imputation_replicate_filename, sep='\t', compression='gzip')

        # subset antigen MM and amino acid MM columns from data frame

        SRTR_agMM = SRTR[["PX_ID", "CAN_RACE", "REC_A_MM_EQUIV_CUR", "REC_B_MM_EQUIV_CUR",
                          "REC_DR_MM_EQUIV_CUR"]]  # "REC_DQ_MM_EQUIV_CUR" is missing from new files

        SRTR_aaMM = SRTR[aamm_in_header]  # Only the AAMM
        # Append the aaMM to the agMM
        SRTR = pd.concat([SRTR_agMM, SRTR_aaMM], axis=1)

        # Convert AgMM to integers and handle NAs due to missing DQ typing
        # 99 - Not tested per SRTR data dictionary
        SRTR['REC_DR_MM_EQUIV_CUR'] = SRTR['REC_DR_MM_EQUIV_CUR'].fillna(99).astype(int)
        SRTR['REC_DQ_MM_EQUIV_CUR'] = SRTR['REC_DQ_MM_EQUIV_CUR'].fillna(99).astype(int)
        SRTR = SRTR.astype({'REC_DR_MM_EQUIV_CUR': 'int', 'REC_DQ_MM_EQUIV_CUR': 'int'})

        if pop == "OVERALL":
            SRTR_pop = SRTR
        else:
            SRTR_pop = SRTR[(SRTR['CAN_RACE'] == pop)]

        # Create high and low FIBERS risk categories with AA positions found
        # Initialize low risk condition, all AA positions are not MM
        low_risk = (SRTR_pop[aamm_in_header[0]] == 0)
        for aamm in aamm_in_header[1:]:  # This allows you to loop through all the amino acids, no matter how many
            low_risk &= (SRTR_pop[aamm] == 0)

        SRTR_Low_Risk = SRTR_pop[low_risk]

        # Initialize high risk condition, at least one of the AA positions are MM
        high_risk = (SRTR_pop[aamm_in_header[0]] >= 1)
        for aamm in aamm_in_header[1:]:
            high_risk |= (SRTR_pop[aamm] >= 1)

        SRTR_High_Risk = SRTR_pop[high_risk]

        SRTR_Overall_Risk = SRTR_pop  # For OVERALL df

        # overall count for each df
        kidney_tx_pair_count = len(SRTR_Overall_Risk)
        FIBERS_Low_Risk_count = len(SRTR_Low_Risk)
        FIBERS_High_Risk_count = len(SRTR_High_Risk)

        # Call the function for each dataframe
        overall_antigen = antigen_count(SRTR_Overall_Risk, replicate, kidney_tx_pair_count)
        OVERALL["Replicate_" + str(replicate)] = overall_antigen
        low_antigen = antigen_count(SRTR_Low_Risk, replicate, FIBERS_Low_Risk_count)
        FIBERS_LOW["Replicate_" + str(replicate)] = low_antigen
        high_antigen = antigen_count(SRTR_High_Risk, replicate, FIBERS_High_Risk_count)
        FIBERS_HIGH["Replicate_" + str(replicate)] = high_antigen

    with pd.ExcelWriter("FIBERS_AgMM_" + pop + "_summary_table.xlsx") as writer:
        OVERALL.to_excel(writer, sheet_name="OVERALL")
        FIBERS_HIGH.to_excel(writer, sheet_name="FIBERS_HIGH")
        FIBERS_LOW.to_excel(writer, sheet_name="FIBERS_LOW")
        
        
# Formatting the Excel sheets for readability
for pop in pops:
    wb = load_workbook(filename="FIBERS_AgMM_" + pop + "_summary_table.xlsx")
    overall_ws = wb['OVERALL']
    high_ws = wb['FIBERS_HIGH']
    low_ws = wb['FIBERS_LOW']

    overall_ws['A1'].value = "Mismatch Category"
    overall_ws['A1'].font = Font(bold=True)
    high_ws['A1'].value = "Mismatch Category"
    high_ws['A1'].font = Font(bold=True)
    low_ws['A1'].value = "Mismatch Category"
    low_ws['A1'].font = Font(bold=True)

    max_column = overall_ws.max_column

    # Increase the width and left-alignment each column
    for column_index in range(1, max_column + 1):
        column_letter = overall_ws.cell(row=1, column=column_index).column_letter
        overall_ws.column_dimensions[column_letter].width = 14
        high_ws.column_dimensions[column_letter].width = 14
        low_ws.column_dimensions[column_letter].width = 14

        for cell in overall_ws[column_letter]:
            cell.alignment = Alignment(horizontal='left')
        for cell in high_ws[column_letter]:
            cell.alignment = Alignment(horizontal='left')
        for cell in low_ws[column_letter]:
            cell.alignment = Alignment(horizontal='left')

    overall_ws.column_dimensions['A'].width = 30
    high_ws.column_dimensions['A'].width = 30
    low_ws.column_dimensions['A'].width = 30

    wb.save(filename="FIBERS_AgMM_" + pop + "_summary_table.xlsx")
