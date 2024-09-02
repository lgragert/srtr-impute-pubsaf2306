###############################################################################
#   SCRIPT NAME:    EpReg_eplet_matrix_generation_array10.py
#   DESCRIPTION:    Generate eplet matrices with all eplet mismatch categories
#   INPUT:          ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_loc.csv
#   OUTPUT:         ./SRTR_output_transfer/eplet_mismatch_Eplet_Matrix_{ARRAY_ID}.csv
#   NOTES:
###############################################################################

import os
import pandas as pd
import ast
import re

# Get the array task ID from the environment variable
ARRAY_ID = int(os.environ['SLURM_ARRAY_TASK_ID'])

# File paths
input_file_path = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_loc.csv'
output_file_path = f'./SRTR_output_transfer/eplet_mismatch_Eplet_Matrix_{ARRAY_ID}.csv'

# Load the dataset
eplet_locus = pd.read_csv(input_file_path)



# Function to extract unique eplets and sort them numerically
def extract_and_sort_eplets(column):
    eplets = set()
    for details in column.dropna():
        try:
            eplets.update(ast.literal_eval(details)) # Convert string to list and update set of eplets
        except:
            pass
    # Sort eplets numerically, considering any preceding letters
    valid_eplets = [eplet for eplet in eplets if re.findall(r'\d+', eplet)]
    eplets = sorted(valid_eplets, key=lambda x: int(re.findall(r'\d+', x)[0]))
    return eplets

# Extract and sort eplets for each locus
unique_eplets_A = extract_and_sort_eplets(eplet_locus['ABC_explanation_A_details'])
unique_eplets_B = extract_and_sort_eplets(eplet_locus['ABC_explanation_B_details'])
unique_eplets_C = extract_and_sort_eplets(eplet_locus['ABC_explanation_C_details'])
unique_eplets_DRB1 = extract_and_sort_eplets(eplet_locus['DRB_explanation_DRB1_details'])
unique_eplets_DRB345 = extract_and_sort_eplets(eplet_locus['DRB_explanation_DRB345_details'])
unique_eplets_DQA1 = extract_and_sort_eplets(eplet_locus['DQ_explanation_DQA1_details'])
unique_eplets_DQB1 = extract_and_sort_eplets(eplet_locus['DQ_explanation_DQB1_details'])

# Create columns for each unique eplet and count occurrences
def add_eplet_columns(df, column, eplets, prefix):
    # Create a dictionary to hold the new columns
    new_columns = {}
    for eplet in eplets:
        # Create a column for each eplet with exact match counts
        # Compute the counts and store them in the dictionary
        new_columns[f'{prefix}_{eplet}'] = df[column].apply(lambda x: ast.literal_eval(x).count(eplet) if pd.notnull(x) else 0)
    # Create a temporary DataFrame from the dictionary (with the same index as the original DataFrame)
    temp_df = pd.DataFrame(new_columns, index=df.index)
    # Concatenate the temporary DataFrame with the original DataFrame
    return pd.concat([df, temp_df], axis=1)

# Add eplet columns to the dataframe
eplet_locus = add_eplet_columns(eplet_locus, 'ABC_explanation_A_details', unique_eplets_A, 'A')
eplet_locus = add_eplet_columns(eplet_locus, 'ABC_explanation_B_details', unique_eplets_B, 'B')
eplet_locus = add_eplet_columns(eplet_locus, 'ABC_explanation_C_details', unique_eplets_C, 'C')
eplet_locus = add_eplet_columns(eplet_locus, 'DRB_explanation_DRB1_details', unique_eplets_DRB1, 'DRB1')
eplet_locus = add_eplet_columns(eplet_locus, 'DRB_explanation_DRB345_details', unique_eplets_DRB345, 'DRB345')
eplet_locus = add_eplet_columns(eplet_locus, 'DQ_explanation_DQA1_details', unique_eplets_DQA1, 'DQA1')
eplet_locus = add_eplet_columns(eplet_locus, 'DQ_explanation_DQB1_details', unique_eplets_DQB1, 'DQB1')



# Drop columns
columns_to_drop = ['ALL_details', 'ALL_explanation',
                   'ABC_quantity', 'ABC_details', 
                   'ABC_explanation_A_details', 'ABC_explanation_A_explanation',
                   'ABC_explanation_B_details', 'ABC_explanation_B_explanation',
                   'ABC_explanation_C_details', 'ABC_explanation_C_explanation',
                   'DRB_quantity', 'DRB_details',
                   'DRB_explanation_DRB1_details', 'DRB_explanation_DRB1_explanation',
                   'DRB_explanation_DRB345_details', 'DRB_explanation_DRB345_explanation',
                   'DQ_quantity', 'DQ_details',
                   'DQ_explanation_DQA1_details', 'DQ_explanation_DQA1_explanation',
                   'DQ_explanation_DQB1_details', 'DQ_explanation_DQB1_explanation',
                   'DP_quantity', 'DP_details', 'DP_explanation',
                   'MICA_quantity', 'MICA_details', 'MICA_explanation']
eplet_matrix = eplet_locus.drop(columns=columns_to_drop, errors='ignore')

# Rename remaining columns
new_column_names = {
    'ALL_quantity': 'total_epmm_count',
    'ABC_explanation_A_quantity': 'A_epmm_count',
    'ABC_explanation_B_quantity': 'B_epmm_count',
    'ABC_explanation_C_quantity': 'C_epmm_count',
    'DRB_explanation_DRB1_quantity': 'DRB1_epmm_count',
    'DRB_explanation_DRB345_quantity': 'DRB345_epmm_count',
    'DQ_explanation_DQA1_quantity': 'DQA1_epmm_count',
    'DQ_explanation_DQB1_quantity': 'DQB1_epmm_count',
}
eplet_matrix = eplet_matrix.rename(columns=new_column_names)

# Save as CSV
eplet_matrix.to_csv(output_file_path, index=False)

print(f"Eplet Matrix saved to {output_file_path}.")
