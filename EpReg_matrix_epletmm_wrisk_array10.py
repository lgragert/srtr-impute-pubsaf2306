###############################################################################
#   SCRIPT NAME:    EpReg_matrix_epletmm_wrisk_array10.py
#   DESCRIPTION:    Combines eplet mismatch and wiebe risk categories
#                   Adds newly generated categories back to original SRTR matrix
#   INPUT:          ./SRTR_AA_MM_9loc_matrix_{ARRAY_ID}.txt.gz
#                   ./SRTR_output_transfer/eplet_mismatch_Eplet_Matrix_{ARRAY_ID}.csv
#                   ./SRTR_output_transfer/eplet_mismatch_Wiebe_Risk_{ARRAY_ID}.csv
#   OUTPUT:         ./SRTR_output_transfer/eplet_mismatch_Matrix_WRisk_{ARRAY_ID}.csv
#                   ./SRTR_AA_MM_9loc_eplet_matrix_{ARRAY_ID}.txt.gz
#   NOTES:
###############################################################################

import os
import pandas as pd

# Get the array task ID from the environment variable
ARRAY_ID = int(os.environ['SLURM_ARRAY_TASK_ID'])

# File paths
input_file_path_SRTR_Matrix = f'./SRTR_AA_MM_9loc_matrix_{ARRAY_ID}.txt.gz'
input_file_path_Eplet_Matrix = f'./SRTR_output_transfer/eplet_mismatch_Eplet_Matrix_{ARRAY_ID}.csv'
input_file_path_Wiebe_Risk = f'./SRTR_output_transfer/eplet_mismatch_Wiebe_Risk_{ARRAY_ID}.csv'
output_file_path_Matrix_WRisk = f'./SRTR_output_transfer/eplet_mismatch_Matrix_WRisk_{ARRAY_ID}.csv'
output_file_path_SRTR_Eplet = f'./SRTR_AA_MM_9loc_eplet_matrix_{ARRAY_ID}.txt.gz'



# Load the datasets
eplet_matrix = pd.read_csv(input_file_path_Eplet_Matrix)
wiebe_risk = pd.read_csv(input_file_path_Wiebe_Risk)

# Select columns to merge from 'Wiebe_Risk'
columns_to_merge = ['PX_ID', 'Max_DR', 'Max_DQ', 'Wiebe_Low', 'Wiebe_Medium', 'Wiebe_High', 'Wiebe_Risk']
wiebe_risk_subset = wiebe_risk[columns_to_merge]

# Merge DataFrames on 'PX_ID'
merged_epletmatrix_wieberisk = pd.merge(wiebe_risk_subset, eplet_matrix, on='PX_ID', how='outer')

# Save as CSV
merged_epletmatrix_wieberisk.to_csv(output_file_path_Matrix_WRisk, index=False)
print(f"Saved eplet matrix (with risk) to {output_file_path_Matrix_WRisk}")



# Read in the original SRTR matrix file 
SRTR_matrix = pd.read_csv(input_file_path_SRTR_Matrix, sep='\t', compression='gzip', dtype={'PX_ID': int})
print(len(SRTR_matrix))
print(SRTR_matrix['PX_ID'].dtype)

# Read in the merged file to handle 'PX_ID' values
eplet_wrisk = pd.read_csv(output_file_path_Matrix_WRisk, dtype={'PX_ID': int})
print(len(eplet_wrisk))
print(eplet_wrisk['PX_ID'].dtype)

# Merge the two dataframes on 'PX_ID'
SRTR_eplet = pd.merge(SRTR_matrix, eplet_wrisk, on='PX_ID', how='left')
print(len(SRTR_eplet))
print(SRTR_eplet['PX_ID'].dtype)

# Save the merged dataframe to a new version of the file (updated with eplet mm categories)
SRTR_eplet.to_csv(output_file_path_SRTR_Eplet, sep='\t', compression='gzip', index=False)
print(f"Saved SRTR eplet matrix to {output_file_path_SRTR_Eplet}")
