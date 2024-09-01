###############################################################################
#   SCRIPT NAME:    EpReg_DRDQ_single_molecule_array10.py
#   DESCRIPTION:    Generate risk categories based on the Wiebe formula
#   INPUT:          ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED.csv
#   OUTPUT:         ./SRTR_output_transfer/eplet_mismatch_Wiebe_Risk_{ARRAY_ID}.csv
#                   ./SRTR_output_transfer/eplet_mismatch_SS_Risk_{ARRAY_ID}.csv
#   NOTES:          Risk logic script modified from hlaR_DRDQ_single_molecule_post_I10.R
###############################################################################

import os
import pandas as pd

# Get the array task ID from the environment variable
ARRAY_ID = int(os.environ['SLURM_ARRAY_TASK_ID'])

# File paths
input_file_path = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED.csv'
output_file_path = f'./SRTR_output_transfer/eplet_mismatch_Wiebe_Risk_{ARRAY_ID}.csv'
output_file_path_ss = f'./SRTR_output_transfer/eplet_mismatch_SS_Risk_{ARRAY_ID}.csv'

# Load in the dataset
eplet_mismatches = pd.read_csv(input_file_path)



# Define risk stratification logic
# High/Med/Low risk from Wiebe et al. AJT 2019
# Figure 3.C values found to be correlated with DSAs and survival

# Calculate max DR mismatch
eplet_mismatches['Max_DR'] = eplet_mismatches[['DRB_explanation_DRB1_allele_1_quantity',
                                               'DRB_explanation_DRB1_allele_2_quantity',
                                               'DRB_explanation_DRB345_allele_1_quantity',
                                               'DRB_explanation_DRB345_allele_2_quantity'
                                               ]].max(axis=1)

# Calculate max DQ mismatch from allele pairs
eplet_mismatches['Max_DQ'] = eplet_mismatches[['DQ_explanation_DQA1_allele_1_quantity',
                                               'DQ_explanation_DQB1_allele_1_quantity'
                                               ]].sum(axis=1).combine(
                                                   eplet_mismatches[['DQ_explanation_DQA1_allele_2_quantity',
                                                                     'DQ_explanation_DQB1_allele_2_quantity'
                                                                     ]].sum(axis=1), max)

# Initialize the risk category columns
eplet_mismatches['Wiebe_Low'] = 0
eplet_mismatches['Wiebe_Medium'] = 0
eplet_mismatches['Wiebe_High'] = 0

# Vectorized risk categorization logic
eplet_mismatches['Wiebe_Low'] = (((eplet_mismatches['Max_DR'] == 0) & (eplet_mismatches['Max_DQ'] == 0))).astype(int)
eplet_mismatches['Wiebe_Medium'] = (((eplet_mismatches['Max_DR'] == 0) & (1 <= eplet_mismatches['Max_DQ']) & (eplet_mismatches['Max_DQ'] <= 8)) |
                                    ((1 <= eplet_mismatches['Max_DR']) & (eplet_mismatches['Max_DR'] <= 6) & (eplet_mismatches['Max_DQ'] <= 8))).astype(int)
eplet_mismatches['Wiebe_High'] = (((eplet_mismatches['Max_DR'] == 0) & (eplet_mismatches['Max_DQ'] >= 9)) |
                                  ((1 <= eplet_mismatches['Max_DR']) & (eplet_mismatches['Max_DR'] <= 6) & (eplet_mismatches['Max_DQ'] >= 9)) |
                                  (eplet_mismatches['Max_DR'] >= 7)).astype(int)

# Generate the risk_category column
eplet_mismatches['Wiebe_Risk'] = pd.Series(['low'] * len(eplet_mismatches))
eplet_mismatches.loc[eplet_mismatches['Wiebe_Medium'] == 1, 'Wiebe_Risk'] = 'interm'
eplet_mismatches.loc[eplet_mismatches['Wiebe_High'] == 1, 'Wiebe_Risk'] = 'high'



# Save the results to a new CSV file
eplet_mismatches.to_csv(output_file_path, index=False)

print(f"Risk categories saved to {output_file_path}.")



# Summary statistics for low/med/high risk categories

# Count the number of occurrences of each risk category
risk_counts = eplet_mismatches['Wiebe_Risk'].value_counts()
print(risk_counts)
risk_fractions = eplet_mismatches['Wiebe_Risk'].value_counts(normalize=True)
print(risk_fractions)

# Combine counts and fractions into a single DataFrame
risk_ss = pd.DataFrame({
    'Counts': risk_counts,
    'Fraction': risk_fractions
})

# Categories (low/med/high) are the index (IOW, part of the df)
# Making the index a column
risk_ss = risk_ss.reset_index()
risk_ss = risk_ss.rename(columns={'index': 'Wiebe_Risk'})

# Save the combined results to a CSV file
risk_ss.to_csv(output_file_path_ss, index=False, header=True)

print(f"Risk summary statistics saved to {output_file_path_ss}.")
