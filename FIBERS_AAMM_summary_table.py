import pandas as pd
from openpyxl.styles import Alignment, Font
from openpyxl import load_workbook
from collections import defaultdict

# Script to generate summary stats for amino acid MM in replicate 1

# Overall Antigen Mismatch Counts Function
def antigen_count(low_high_bin, total_count, high_or_low, overall):
    # Initialize every row
    antigen_DR_NA_count = len(low_high_bin[low_high_bin['REC_DR_MM_EQUIV_CUR'] == 99])
    antigen_DQ_NA_count = len(low_high_bin[low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 99])

    antigen_ABDR_0MM_count = len(low_high_bin[(low_high_bin['REC_A_MM_EQUIV_CUR'] == 0) &
                                              (low_high_bin['REC_B_MM_EQUIV_CUR'] == 0) &
                                              (low_high_bin['REC_DR_MM_EQUIV_CUR'] == 0)])

    antigen_DR_0MM_count = len(low_high_bin[low_high_bin['REC_DR_MM_EQUIV_CUR'] == 0])
    antigen_DQ_0MM_count = len(low_high_bin[low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 0])

    antigen_DR_1MM_count = len(low_high_bin[low_high_bin['REC_DR_MM_EQUIV_CUR'] == 1])
    antigen_DQ_1MM_count = len(low_high_bin[low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 1])

    antigen_DR_2MM_count = len(low_high_bin[low_high_bin['REC_DR_MM_EQUIV_CUR'] == 2])
    antigen_DQ_2MM_count = len(low_high_bin[low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 2])

    antigen_DR_0MM_DQ_0MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 0)])

    antigen_DR_0MM_DQ_1MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 1)])

    antigen_DR_1MM_DQ_0MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 0)])

    antigen_DR_1MM_DQ_1MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 1)])

    antigen_DR_0MM_DQ_2MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 2)])

    antigen_DR_2MM_DQ_0MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 0)])

    antigen_DR_1MM_DQ_2MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 2)])

    antigen_DR_2MM_DQ_1MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 1)])

    antigen_DR_2MM_DQ_2MM_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                   (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 2)])

    antigen_DR_0MM_DQ_NA_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 0) &
                                                  (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 99)])

    antigen_DR_1MM_DQ_NA_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 1) &
                                                  (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 99)])

    antigen_DR_2MM_DQ_NA_count = len(low_high_bin[(low_high_bin['REC_DR_MM_EQUIV_CUR'] == 2) &
                                                  (low_high_bin['REC_DQ_MM_EQUIV_CUR'] == 99)])
    # Make the rows depending on if overall count or not
    if overall == True:
        antigen_DR_NA_perc = str(antigen_DR_NA_count) + " (" + str(format(100 * antigen_DR_NA_count / total_count, '.2f')) + "%)"
        antigen_DQ_NA_perc = str(antigen_DQ_NA_count) + " (" + str(format(100 * antigen_DQ_NA_count / total_count, ".2f")) + "%)"
        antigen_ABDR_0MM_perc = str(antigen_ABDR_0MM_count) + " (" + str(format(100 * antigen_ABDR_0MM_count / total_count, '.2f')) + "%)"
        antigen_DR_0MM_perc = str(antigen_DR_0MM_count) + " (" + str(format(100 * antigen_DR_0MM_count / total_count, '.2f')) + "%)"
        antigen_DQ_0MM_perc = str(antigen_DQ_0MM_count) + " (" + str(format(100 * antigen_DQ_0MM_count / total_count, '.2f')) + "%)"
        antigen_DR_1MM_perc = str(antigen_DR_1MM_count) + " (" + str(format(100 * antigen_DR_1MM_count / total_count, '.2f')) + "%)"
        antigen_DQ_1MM_perc = str(antigen_DQ_1MM_count) + " (" + str(format(100 * antigen_DQ_1MM_count / total_count, '.2f')) + "%)"
        antigen_DR_2MM_perc = str(antigen_DR_2MM_count) + " (" + str(format(100 * antigen_DR_2MM_count / total_count, '.2f')) + "%)"
        antigen_DQ_2MM_perc = str(antigen_DQ_2MM_count) + " (" + str(format(100 * antigen_DQ_2MM_count / total_count, '.2f')) + "%)"
        antigen_DR_0MM_DQ_0MM_perc = str(antigen_DR_0MM_DQ_0MM_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_0MM_count / total_count, '.2f')) + "%)"
        antigen_DR_0MM_DQ_1MM_perc = str(antigen_DR_0MM_DQ_1MM_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_1MM_count / total_count, '.2f')) + "%)"
        antigen_DR_1MM_DQ_0MM_perc = str(antigen_DR_1MM_DQ_0MM_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_0MM_count / total_count, '.2f')) + "%)"
        antigen_DR_1MM_DQ_1MM_perc = str(antigen_DR_1MM_DQ_1MM_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_1MM_count / total_count, '.2f')) + "%)"
        antigen_DR_0MM_DQ_2MM_perc = str(antigen_DR_0MM_DQ_2MM_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_2MM_count / total_count, '.2f')) + "%)"
        antigen_DR_2MM_DQ_0MM_perc = str(antigen_DR_2MM_DQ_0MM_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_0MM_count / total_count, '.2f')) + "%)"
        antigen_DR_1MM_DQ_2MM_perc = str(antigen_DR_1MM_DQ_2MM_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_2MM_count / total_count, '.2f')) + "%)"
        antigen_DR_2MM_DQ_1MM_perc = str(antigen_DR_2MM_DQ_1MM_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_1MM_count / total_count, '.2f')) + "%)"
        antigen_DR_2MM_DQ_2MM_perc = str(antigen_DR_2MM_DQ_2MM_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_2MM_count / total_count, '.2f')) + "%)"
        antigen_DR_0MM_DQ_NA_perc = str(antigen_DR_0MM_DQ_NA_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_NA_count / total_count, '.2f')) + "%)"
        antigen_DR_1MM_DQ_NA_perc = str(antigen_DR_1MM_DQ_NA_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_NA_count / total_count, '.2f')) + "%)"
        antigen_DR_2MM_DQ_NA_perc = str(antigen_DR_2MM_DQ_NA_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_NA_count / total_count, '.2f')) + "%)"
    else:
        total_count[0].replace(to_replace="0 ", value="1", inplace=True)  # Cannot divide by zero, so replace with one
        antigen_DR_NA_perc = str(antigen_DR_NA_count) + " (" + str(format(100 * antigen_DR_NA_count / int(total_count.loc['Missing-DR AgMM Count'][0]), '.2f')) + "%)"
        antigen_DQ_NA_perc = str(antigen_DQ_NA_count) + " (" + str(format(100 * antigen_DQ_NA_count / int(total_count.loc['Missing-DQ AgMM Count'][0]), ".2f")) + "%)"
        antigen_ABDR_0MM_perc = str(antigen_ABDR_0MM_count) + " (" + str(format(100 * antigen_ABDR_0MM_count / int(total_count.loc['0-ABDR AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_0MM_perc = str(antigen_DR_0MM_count) + " (" + str(format(100 * antigen_DR_0MM_count / int(total_count.loc['0-DR AgMM Count'][0]), '.2f')) + "%)"
        antigen_DQ_0MM_perc = str(antigen_DQ_0MM_count) + " (" + str(format(100 * antigen_DQ_0MM_count / int(total_count.loc['0-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_1MM_perc = str(antigen_DR_1MM_count) + " (" + str(format(100 * antigen_DR_1MM_count / int(total_count.loc['1-DR AgMM Count'][0]), '.2f')) + "%)"
        antigen_DQ_1MM_perc = str(antigen_DQ_1MM_count) + " (" + str(format(100 * antigen_DQ_1MM_count / int(total_count.loc['1-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_2MM_perc = str(antigen_DR_2MM_count) + " (" + str(format(100 * antigen_DR_2MM_count / int(total_count.loc['2-DR AgMM Count'][0]), '.2f')) + "%)"
        antigen_DQ_2MM_perc = str(antigen_DQ_2MM_count) + " (" + str(format(100 * antigen_DQ_2MM_count / int(total_count.loc['2-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_0MM_DQ_0MM_perc = str(antigen_DR_0MM_DQ_0MM_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_0MM_count / int(total_count.loc['0-DR and 0-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_0MM_DQ_1MM_perc = str(antigen_DR_0MM_DQ_1MM_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_1MM_count / int(total_count.loc['0-DR and 1-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_1MM_DQ_0MM_perc = str(antigen_DR_1MM_DQ_0MM_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_0MM_count / int(total_count.loc['1-DR and 0-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_1MM_DQ_1MM_perc = str(antigen_DR_1MM_DQ_1MM_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_1MM_count / int(total_count.loc['1-DR and 1-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_0MM_DQ_2MM_perc = str(antigen_DR_0MM_DQ_2MM_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_2MM_count / int(total_count.loc['0-DR and 2-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_2MM_DQ_0MM_perc = str(antigen_DR_2MM_DQ_0MM_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_0MM_count / int(total_count.loc['2-DR and 0-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_1MM_DQ_2MM_perc = str(antigen_DR_1MM_DQ_2MM_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_2MM_count / int(total_count.loc['1-DR and 2-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_2MM_DQ_1MM_perc = str(antigen_DR_2MM_DQ_1MM_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_1MM_count / int(total_count.loc['2-DR and 1-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_2MM_DQ_2MM_perc = str(antigen_DR_2MM_DQ_2MM_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_2MM_count / int(total_count.loc['2-DR and 2-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_0MM_DQ_NA_perc = str(antigen_DR_0MM_DQ_NA_count) + " (" + str(format(100 * antigen_DR_0MM_DQ_NA_count / int(total_count.loc['0-DR and Missing-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_1MM_DQ_NA_perc = str(antigen_DR_1MM_DQ_NA_count) + " (" + str(format(100 * antigen_DR_1MM_DQ_NA_count / int(total_count.loc['1-DR and Missing-DQ AgMM Count'][0]), '.2f')) + "%)"
        antigen_DR_2MM_DQ_NA_perc = str(antigen_DR_2MM_DQ_NA_count) + " (" + str(format(100 * antigen_DR_2MM_DQ_NA_count / int(total_count.loc['2-DR and Missing-DQ AgMM Count'][0]), '.2f')) + "%)"

    if high_or_low == "high":
        antigen_list = [" ", " ", "", antigen_DR_NA_perc, antigen_DQ_NA_perc, antigen_ABDR_0MM_perc, antigen_DR_0MM_perc,
                        antigen_DQ_0MM_perc, antigen_DR_1MM_perc, antigen_DQ_1MM_perc, antigen_DR_2MM_perc,
                        antigen_DQ_2MM_perc, antigen_DR_0MM_DQ_0MM_perc, antigen_DR_0MM_DQ_1MM_perc,
                        antigen_DR_1MM_DQ_0MM_perc, antigen_DR_1MM_DQ_1MM_perc,
                        antigen_DR_0MM_DQ_2MM_perc, antigen_DR_2MM_DQ_0MM_perc,
                        antigen_DR_1MM_DQ_2MM_perc, antigen_DR_2MM_DQ_1MM_perc,
                        antigen_DR_2MM_DQ_2MM_perc, antigen_DR_0MM_DQ_NA_perc,
                        antigen_DR_1MM_DQ_NA_perc, antigen_DR_2MM_DQ_NA_perc]
        # "Leave blank"'s are due to formatting issues, but they will be fixed later in the process
        antigen_key = ["Leave blank", "leave blank", "AA-MM within Ag-MM and FIBERS High Risk Combos",
                       "Missing-DR AgMM Count", "Missing-DQ AgMM Count", "0-ABDR AgMM Count", "0-DR AgMM Count",
                       "0-DQ AgMM Count", "1-DR AgMM Count", "1-DQ AgMM Count", "2-DR AgMM Count", "2-DQ AgMM Count",
                       "0-DR and 0-DQ AgMM Count", "0-DR and 1-DQ AgMM Count", "1-DR and 0-DQ AgMM Count",
                       "1-DR and 1-DQ AgMM Count", "0-DR and 2-DQ AgMM Count", "2-DR and 0-DQ AgMM Count",
                       "1-DR and 2-DQ AgMM Count", "2-DR and 1-DQ AgMM Count", "2-DR and 2-DQ AgMM Count",
                       "0-DR and Missing-DQ AgMM Count", "1-DR and Missing-DQ AgMM Count", "2-DR and Missing-DQ AgMM Count"]
    if high_or_low == "low":
        antigen_list = ["", antigen_DR_NA_perc, antigen_DQ_NA_perc, antigen_ABDR_0MM_perc,
                        antigen_DR_0MM_perc, antigen_DQ_0MM_perc, antigen_DR_1MM_perc,
                        antigen_DQ_1MM_perc, antigen_DR_2MM_perc, antigen_DQ_2MM_perc,
                        antigen_DR_0MM_DQ_0MM_perc, antigen_DR_0MM_DQ_1MM_perc,
                        antigen_DR_1MM_DQ_0MM_perc, antigen_DR_1MM_DQ_1MM_perc,
                        antigen_DR_0MM_DQ_2MM_perc, antigen_DR_2MM_DQ_0MM_perc,
                        antigen_DR_1MM_DQ_2MM_perc, antigen_DR_2MM_DQ_1MM_perc,
                        antigen_DR_2MM_DQ_2MM_perc, antigen_DR_0MM_DQ_NA_perc,
                        antigen_DR_1MM_DQ_NA_perc, antigen_DR_2MM_DQ_NA_perc]
        antigen_key = ["AA-MM within Ag-MM and FIBERS Low Risk Combos", "Missing-DR AgMM Count",
                       "Missing-DQ AgMM Count", "0-ABDR AgMM Count", "0-DR AgMM Count", "0-DQ AgMM Count",
                       "1-DR AgMM Count", "1-DQ AgMM Count", "2-DR AgMM Count", "2-DQ AgMM Count",
                       "0-DR and 0-DQ AgMM Count", "0-DR and 1-DQ AgMM Count", "1-DR and 0-DQ AgMM Count",
                       "1-DR and 1-DQ AgMM Count", "0-DR and 2-DQ AgMM Count", "2-DR and 0-DQ AgMM Count",
                       "1-DR and 2-DQ AgMM Count", "2-DR and 1-DQ AgMM Count", "2-DR and 2-DQ AgMM Count",
                       "0-DR and Missing-DQ AgMM Count", "1-DR and Missing-DQ AgMM Count", "2-DR and Missing-DQ AgMM Count"]

    antigen_dictionary = dict(zip(antigen_key, antigen_list))

    return antigen_dictionary


# Specific Amino Acids Count Function: DRB1_13 _26 and DQB1_30 _55
def specific_amino_acid(SRTR, high_list, low_list, high_overall, low_overall, amino_acid):

    # Initialize SRTR to have specific AA-MM present
    col_AA = 'MM_' + amino_acid
    SRTR_AA = SRTR[(SRTR[col_AA] >= 1)]

    high_aa_ag_dict = defaultdict(lambda: defaultdict(dict))
    low_aa_ag_dict = defaultdict(lambda: defaultdict(dict))
    total_low_count_dict = defaultdict()
    total_high_count_dict = defaultdict()

    for idx, high in enumerate(high_list):
        # Checks to compare AA-MM presence in bin with high risk
        high_bin_SRTR = high[high[col_AA].isin(SRTR_AA[col_AA])]
        total_count = len(high)
        high_count = len(high_bin_SRTR)
        binSRTR_count = str(high_count) + " (" + str(format(100 * high_count / total_count, ".2f")) + "%)"
        total_high_count_dict["FIBERS_Bin_" + str(idx + 1)] = high_count
        # Need these totals for the denominators for the Ag-MM Counts
        no_percents_high = high_overall["FIBERS_Bin_" + str(idx + 1)].str.split("(", expand=True)
        no_percents_high = no_percents_high.drop(columns=1, index=["Leave blank", "leave blank",
                                                                   "AA-MM within Ag-MM and FIBERS High Risk Combos",
                                                                   "FIBERS High Risk >=1-AA-MM"])
        antigen_dict = antigen_count(high_bin_SRTR, no_percents_high, "high", False)
        high_aa_ag_dict["FIBERS_Bin_" + str(idx + 1)] = antigen_dict
        high_aa_ag_dict["FIBERS_Bin_" + str(idx + 1)]['FIBERS High Risk >=1-AA-MM'] = binSRTR_count

    for idx, low in enumerate(low_list):
        # Checks to compare AA-MM presence in bin with high risk
        low_bin_SRTR = low[low[col_AA].isin(SRTR_AA[col_AA])]
        total_count = len(low)
        low_count = len(low_bin_SRTR)
        binSRTR_count = str(low_count) + " (" + str(format(100 * low_count / total_count, ".2f")) + "%)"
        total_low_count_dict["FIBERS_Bin_" + str(idx + 1)] = low_count
        # Need these totals for the denominators for the Ag-MM Counts
        no_percents_low = low_overall["FIBERS_Bin_" + str(idx + 1)].str.split("(", expand=True)
        no_percents_low = no_percents_low.drop(columns=1, index=["AA-MM within Ag-MM and FIBERS Low Risk Combos",
                                                                 "FIBERS Low Risk 0-AA-MM"])
        antigen_dict = antigen_count(low_bin_SRTR, no_percents_low, "low", False)
        low_aa_ag_dict["FIBERS_Bin_" + str(idx + 1)] = antigen_dict
        low_aa_ag_dict["FIBERS_Bin_" + str(idx + 1)]['FIBERS Low Risk 0-AA-MM'] = binSRTR_count

    aa_high = pd.DataFrame.from_dict(high_aa_ag_dict, orient='columns')
    aa_low = pd.DataFrame.from_dict(low_aa_ag_dict, orient='columns')
    # Add a total AA-MM Count at the top of the spreadsheet
    total_low_count_df = pd.DataFrame.from_dict(total_low_count_dict, orient='index', columns=["Low Risk"])
    total_high_count_df = pd.DataFrame.from_dict(total_high_count_dict, orient='index', columns=["High Risk"])
    total_risk_count = pd.concat([total_low_count_df, total_high_count_df], axis=1)
    total_risk_count_df = pd.DataFrame()
    total_risk_count_df['Total_MM'] = total_risk_count.sum(axis=1)
    total_risk_count_df = total_risk_count_df.transpose()

    aa_ag_df = pd.concat([aa_high, aa_low])
    aa_ag_df = pd.concat([total_risk_count_df, aa_ag_df])

    aa_ag_df.name = amino_acid  # Will Create a Nested Dictionary with Named DataFrames

    return aa_ag_df
    
# Reference table of AA-MM in each bin
table_df = pd.DataFrame(data={"DRB1": ["11, 26, 37, 77", "12", "13, 26, 77", "12, 14", "13, 26", "13, 26", "10, 13, 32, 73", "13, 33", "13, 26", "13, 26"], 
                              "DQB1": ["53, 89", "55, 71, 77", "30", "55, 74", "30, 55, 89", "30, 55", "30, 55", "30, 55, 74", "30, 55", "55"],
                              "A": ["", "", "", "17", "", "", "", "", "", ""]}, index=range(1, 11))

AAMM_tabs = ['DRB1_13', 'DRB1_26', 'DQB1_30', 'DQB1_55', 'DRB1_71', 'DRB1_74', 'DRB1_86', 'DRB1_60', 'DRB1_73']

# SRTR using only one replicate/realization
replicate = 1
SRTR_imputation_replicate_filename = "SRTR_AA_MM_matrix_grffail_" + str(replicate) + ".txt"
SRTR_df = pd.read_csv(SRTR_imputation_replicate_filename, sep='\t')

# Subset amino acid MM and antigen MM columns from Dataframe
SRTR_df = SRTR_df[["PX_ID", "CAN_RACE", "REC_A_MM_EQUIV_CUR", "REC_B_MM_EQUIV_CUR", "REC_DR_MM_EQUIV_CUR",
                   "REC_DQ_MM_EQUIV_CUR", "MM_DRB1_13", "MM_DRB1_26", "MM_DQB1_30", "MM_DQB1_55", "MM_DRB1_11",
                   "MM_DRB1_37", "MM_DQB1_53", "MM_DRB1_77", "MM_DQB1_89", 'MM_DRB1_12', 'MM_DQB1_71', 'MM_DQB1_77',
                   'MM_A_17', 'MM_DRB1_14', 'MM_DQB1_74', 'MM_DRB1_10', 'MM_DRB1_32', 'MM_DRB1_73', 'MM_DRB1_33',
                   'MM_DRB1_71', 'MM_DRB1_74', 'MM_DRB1_86', 'MM_DRB1_60']]

# Convert AgMM to integers and handle NAs due to missing DQ typing
# 99 - Not tested per SRTR_df data dictionary
SRTR_df['REC_DR_MM_EQUIV_CUR'] = SRTR_df['REC_DR_MM_EQUIV_CUR'].fillna(99).astype(int)
SRTR_df['REC_DQ_MM_EQUIV_CUR'] = SRTR_df['REC_DQ_MM_EQUIV_CUR'].fillna(99).astype(int)
SRTR_df = SRTR_df.astype({'REC_DR_MM_EQUIV_CUR': 'int', 'REC_DQ_MM_EQUIV_CUR': 'int'})


# Loop through each population and overall
pops = ["OVERALL", "AFA", "CAU", "HIS", "HPI", "NAM", "ASI", "Multi-Racial"]

# Initialize all the bins and make nested dictionaries for each sheet
high_count_dict = defaultdict(lambda: defaultdict(dict))
low_count_dict = defaultdict(lambda: defaultdict(dict))
amino_acid_dict = defaultdict(lambda: defaultdict(dict))
AAMM_df = defaultdict(lambda: defaultdict(dict))

for pop in pops:
    print("Working on population: ", pop)

    if pop == "OVERALL":
        SRTR = SRTR_df
    else:
        SRTR = SRTR_df[(SRTR_df['CAN_RACE'] == pop)]

    overall_SRTR_count = len(SRTR)
    print("Overall Count:", overall_SRTR_count)

    lowfreq_bin_1 = SRTR[(SRTR['MM_DRB1_11'] == 0) & (SRTR['MM_DRB1_26'] == 0) &
                         (SRTR['MM_DRB1_37'] == 0) & (SRTR['MM_DRB1_77'] == 0) &
                         (SRTR['MM_DQB1_53'] == 0) & (SRTR['MM_DQB1_89'] == 0)]
    highfreq_bin_1 = SRTR[(SRTR['MM_DRB1_11'] >= 1) | (SRTR['MM_DRB1_26'] >= 1) |
                          (SRTR['MM_DRB1_37'] >= 1) | (SRTR['MM_DRB1_77'] >= 1) |
                          (SRTR['MM_DQB1_53'] >= 1) | (SRTR['MM_DQB1_89'] >= 1)]

    lowfreq_bin_2 = SRTR[(SRTR['MM_DRB1_12'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DQB1_71'] == 0) & (SRTR['MM_DQB1_77'] == 0)]
    highfreq_bin_2 = SRTR[(SRTR['MM_DRB1_12'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DQB1_71'] >= 1) | (SRTR['MM_DQB1_77'] >= 1)]

    lowfreq_bin_3 = SRTR[(SRTR['MM_DRB1_13'] == 0) & (SRTR['MM_DRB1_26'] == 0) &
                         (SRTR['MM_DRB1_77'] == 0) & (SRTR['MM_DQB1_30'] == 0)]
    highfreq_bin_3 = SRTR[(SRTR['MM_DRB1_13'] >= 1) | (SRTR['MM_DRB1_26'] >= 1) |
                          (SRTR['MM_DRB1_77'] >= 1) | (SRTR['MM_DQB1_30'] >= 1)]

    lowfreq_bin_4 = SRTR[(SRTR['MM_A_17'] == 0) & (SRTR['MM_DRB1_12'] == 0) &
                         (SRTR['MM_DRB1_14'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DQB1_74'] == 0)]
    highfreq_bin_4 = SRTR[(SRTR['MM_A_17'] >= 1) | (SRTR['MM_DRB1_12'] >= 1) |
                          (SRTR['MM_DRB1_14'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DQB1_74'] >= 1)]

    lowfreq_bin_5 = SRTR[(SRTR['MM_DRB1_13'] == 0) & (SRTR['MM_DRB1_26'] == 0) &
                         (SRTR['MM_DQB1_30'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DQB1_89'] == 0)]
    highfreq_bin_5 = SRTR[(SRTR['MM_DRB1_13'] >= 1) | (SRTR['MM_DRB1_26'] >= 1) |
                          (SRTR['MM_DQB1_30'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DQB1_89'] >= 1)]

    lowfreq_bin_6 = SRTR[(SRTR['MM_DQB1_30'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DRB1_13'] == 0) & (SRTR['MM_DRB1_26'] == 0)]
    highfreq_bin_6 = SRTR[(SRTR['MM_DQB1_30'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DRB1_13'] >= 1) | (SRTR['MM_DRB1_26'] >= 1)]

    lowfreq_bin_7 = SRTR[(SRTR['MM_DQB1_30'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DRB1_10'] == 0) & (SRTR['MM_DRB1_13'] == 0) &
                         (SRTR['MM_DRB1_32'] == 0) & (SRTR['MM_DRB1_73'] == 0)]
    highfreq_bin_7 = SRTR[(SRTR['MM_DQB1_30'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DRB1_10'] >= 1) | (SRTR['MM_DRB1_13'] >= 1) |
                          (SRTR['MM_DRB1_32'] >= 1) | (SRTR['MM_DRB1_73'] >= 1)]

    lowfreq_bin_8 = SRTR[(SRTR['MM_DQB1_30'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DQB1_74'] == 0) & (SRTR['MM_DRB1_13'] == 0) &
                         (SRTR['MM_DRB1_33'] == 0)]
    highfreq_bin_8 = SRTR[(SRTR['MM_DQB1_30'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DQB1_74'] >= 1) | (SRTR['MM_DRB1_13'] >= 1) |
                          (SRTR['MM_DRB1_33'] >= 1)]

    lowfreq_bin_9 = SRTR[(SRTR['MM_DQB1_30'] == 0) & (SRTR['MM_DQB1_55'] == 0) &
                         (SRTR['MM_DRB1_13'] == 0) & (SRTR['MM_DRB1_26'] == 0)]
    highfreq_bin_9 = SRTR[(SRTR['MM_DQB1_30'] >= 1) | (SRTR['MM_DQB1_55'] >= 1) |
                          (SRTR['MM_DRB1_13'] >= 1) | (SRTR['MM_DRB1_26'] >= 1)]

    lowfreq_bin_10 = SRTR[(SRTR['MM_DQB1_55'] == 0) & (SRTR['MM_DRB1_13'] == 0) &
                          (SRTR['MM_DRB1_26'] == 0)]
    highfreq_bin_10 = SRTR[(SRTR['MM_DQB1_55'] >= 1) | (SRTR['MM_DRB1_13'] >= 1) |
                           (SRTR['MM_DRB1_26'] >= 1)]

    high_bins = [highfreq_bin_1, highfreq_bin_2, highfreq_bin_3, highfreq_bin_4, highfreq_bin_5, highfreq_bin_6, highfreq_bin_7, highfreq_bin_8, highfreq_bin_9, highfreq_bin_10]

    low_bins = [lowfreq_bin_1, lowfreq_bin_2, lowfreq_bin_3, lowfreq_bin_4, lowfreq_bin_5, lowfreq_bin_6, lowfreq_bin_7, lowfreq_bin_8, lowfreq_bin_9, lowfreq_bin_10]

    # For OVERALL Excel Sheet
    for idx, high in enumerate(high_bins):
        high_total = len(high)
        high_overall_count = str(high_total) + " (" + str(format(100 * high_total / overall_SRTR_count, ".2f")) + "%)"
        antigen_high = antigen_count(high, high_total, "high", True)
        high_count_dict["FIBERS_Bin_" + str(idx + 1)] = antigen_high
        high_count_dict["FIBERS_Bin_" + str(idx + 1)]["FIBERS High Risk >=1-AA-MM"] = high_overall_count

    for idx, low in enumerate(low_bins):
        low_total = len(low)
        low_overall_count = str(low_total) + " (" + str(format(100 * low_total / overall_SRTR_count, ".2f")) + "%)"
        antigen_low = antigen_count(low, low_total, "low", True)
        low_count_dict["FIBERS_Bin_" + str(idx + 1)] = antigen_low
        low_count_dict["FIBERS_Bin_" + str(idx + 1)]["FIBERS Low Risk 0-AA-MM"] = low_overall_count

    overall_high = pd.DataFrame.from_dict(high_count_dict, orient='columns')
    overall_low = pd.DataFrame.from_dict(low_count_dict, orient='columns')
    OVERALL = pd.concat([overall_high, overall_low])
    # Adds a row for the total population in each bin
    pop_row = pd.DataFrame({"FIBERS_Bin_1": overall_SRTR_count, "FIBERS_Bin_2": overall_SRTR_count, "FIBERS_Bin_3": overall_SRTR_count, "FIBERS_Bin_4": overall_SRTR_count,
                            "FIBERS_Bin_5": overall_SRTR_count, "FIBERS_Bin_6": overall_SRTR_count, "FIBERS_Bin_7": overall_SRTR_count, "FIBERS_Bin_8": overall_SRTR_count,
                            "FIBERS_Bin_9": overall_SRTR_count, "FIBERS_Bin_10": overall_SRTR_count}, index=["Total " + pop + " Population"])
    OVERALL = pd.concat([pop_row, OVERALL])

    # For Each Single AA-MM Excel Sheet
    for AA in AAMM_tabs:
        AAMM_df[AA] = specific_amino_acid(SRTR, high_bins, low_bins, overall_high, overall_low, AA)

    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_colwidth', None)
    # print(AAMM_df)

    with pd.ExcelWriter("FIBERS_AAMM_" + pop + "_summary_table.xlsx") as writer:
        OVERALL.to_excel(writer, sheet_name="OVERALL")
        for AA in AAMM_tabs:
            AAMM_df[AA].to_excel(writer, sheet_name=AA + "_COUNT")
        table_df.to_excel(writer, sheet_name="AA-MM MEMBERSHIP")


# For Excel sheets
FIBERS_heading = NamedStyle(name="FIBERS_heading")
FIBERS_heading.font = Font(bold=True, color="00FF0000")

# Function for what each Excel sheet should look like
def excel_adjuster(worksheet, which_sheet):

    max_column = worksheet.max_column

    for column_index in range(1, max_column + 1):
        column_letter = worksheet.cell(row=1, column=column_index).column_letter
        worksheet.column_dimensions[column_letter].width = 15

        for cell in worksheet[column_letter]:
            cell.alignment = Alignment(horizontal='left')

    if which_sheet == "TABLE":
        worksheet['A1'].value = "Top Bin For Each Imputation"
        worksheet['A1'].style = FIBERS_heading
        worksheet.column_dimensions['A'].width = 25
        for column_index in range(1, max_column + 1):
            column_letter = worksheet.cell(row=1, column=column_index).column_letter
            for cell in worksheet[column_letter]:
                cell.alignment = Alignment(horizontal='center')
    elif which_sheet == "OVERALL":
        worksheet['A1'].value = "MM Category"
    else:
        worksheet['A1'].value = "AA-MM within FIBERS MM Category"

    if which_sheet != "TABLE":
        worksheet.column_dimensions['A'].width = 38
        worksheet.move_range("A27:K27", rows=-24)
        worksheet.move_range("A50:K50", rows=-46)
        worksheet.move_range("A28:K49", rows=-1)
        worksheet['A1'].style = FIBERS_heading
        worksheet['A5'].style = FIBERS_heading
        worksheet['A27'].style = FIBERS_heading

    return worksheet


# Formatting the Excel sheets for readability
for pop in pops:
    wb = load_workbook(filename="FIBERS_AAMM_" + pop + "_summary_table.xlsx")
    wb.add_named_style(FIBERS_heading)
    overall_ws = wb['OVERALL']
    table_ws = wb['AA-MM MEMBERSHIP']

    excel_adjuster(overall_ws, "OVERALL")
    excel_adjuster(table_ws, "TABLE")
    for AA in AAMM_tabs:
        excel_adjuster(wb[AA + "_COUNT"], "AA-MM")

    # Create a row for a new heading for total population
    overall_ws.insert_rows(2)
    overall_ws['A2'].value = "Total Population"
    overall_ws['A2'].style = FIBERS_heading
    if pop == "OVERALL":
        overall_ws['A3'].value = "Total Population in SRTR Kidney"
        overall_ws.font = Font(bold=True)
    else:
        overall_ws['A3'].value = "Total " + pop + " Population"
        overall_ws.font = Font(bold=True)

    for AA in AAMM_tabs:
        aa_ws = wb[AA + "_COUNT"]
        aa_ws.insert_rows(2)
        aa_ws['A2'].value = "Total " + AA + "_AAMM"
        aa_ws['A2'].style = FIBERS_heading
        if pop == "OVERALL":
            aa_ws['A3'].value = "Total Population with " + AA + "_AAMM"
            aa_ws.font = Font(bold=True)
        else:
            aa_ws['A3'].value = "Total " + pop + " Population with " + AA + "_AAMM"
            aa_ws.font = Font(bold=True)

    wb.save(filename="FIBERS_AAMM_" + pop + "_summary_table.xlsx")





