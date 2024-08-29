###############################################################################
#   SCRIPT NAME:    eplet_chunk_processing_array.py
#   DESCRIPTION:    Script to process eplet chunks returned by the Eplet Registry API
#   INPUT:          ./SRTR_output_API_{ARRAY_ID}.tar.gz
#   OUTPUT:         ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_RAW.csv
#                   ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED.csv
#                   ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_labels.csv
#                   ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_loc.csv
#                   ./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_loc_labels.csv
#   NOTES:
###############################################################################

import os
import tarfile
import pandas as pd
import ast



# Get the array task ID from the environment variable
ARRAY_ID = int(os.environ['SLURM_ARRAY_TASK_ID'])

# File paths
input_file_path = f'./SRTR_output_API_{ARRAY_ID}.tar.gz'
output_file_path_RAW = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_RAW.csv'
output_file_path_PROC_ale = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED.csv'
output_file_path_PROC_ale_labels = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_labels.csv'
output_file_path_PROC_loc = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_loc.csv'
output_file_path_PROC_loc_labels = f'./SRTR_output_processed_{ARRAY_ID}/hla_eplet_calculator_output_{ARRAY_ID}_PROCESSED_loc_labels.csv'



####################

# Open the tar.gz file
with tarfile.open('SRTR_output_API_{ARRAY_ID}.tar.gz', 'r:gz') as tar:

    # Only pull/filter out .csv files and
    # Sort the file names numerically based on the chunk number
    file_names = sorted(
        [name for name in tar.getnames() if name.endswith('.csv')],
        key=lambda x: int(x.split('_')[-1].split('.')[0])
        )

    # Initialize an empty list to store DataFrames
    df_list = []

    # Iterate over each sorted file
    for file_name in file_names:
            extracted_file = tar.extractfile(file_name)
            df = pd.read_csv(extracted_file)
            df_list.append(df)

    # Concatenate all DataFrames into one large DataFrame
    eplet_dataframe = pd.concat(df_list, ignore_index=True)

# Save the combined DataFrame as a single CSV file 
eplet_dataframe.to_csv(output_file_path_RAW, index=False, sep=',')

####################

# Reorganize the dataframe

# 'PX_ID' at beginning of dataframe
PXIDcol = eplet_dataframe.pop('PX_ID')
eplet_dataframe.insert(0, 'PX_ID', PXIDcol)

# Group columns by name and re-order them

# Define a grouping order for the columns
group_order = ['PX', 'ALL', 'ABC',
               'DRB_quantity', 'DRB_details',
               'DRB_explanation_DRB1', 'DRB_explanation_DRB345',
               'DQ', 'DP', 'MICA']

# Group and reorder columns
ordered_columns = []
for group in group_order:
    grouped_cols = [col for col in eplet_dataframe.columns if col.startswith(group)]
    ordered_columns.extend(grouped_cols)

# Reindex the DataFrame columns according to the defined order
eplet_dataframe = eplet_dataframe[ordered_columns]

# Processed by allele 1 and allele 2
eplet_dataframe.to_csv(output_file_path_PROC_ale, index=False, sep=',')

# Removing the prefix '{prefix}_explanation_' from column names
eplet_dataframe.columns = eplet_dataframe.columns.str.replace('ABC_explanation_', '').str.replace('DRB_explanation_', '').str.replace('DQ_explanation_', '')

# Cleaned version
eplet_dataframe.to_csv(output_file_path_PROC_ale_labels, index=False, sep=',')

####################

# Split by locus

# Reset eplet dataframe
eplet_dataframe = pd.read_csv(output_file_path_PROC_ale)

# For each locus:
Locus = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1']
ABC_loci = ['A', 'B', 'C']
DR_loci = ['DRB1', 'DRB345']
DQ_loci = ['DQA1', 'DQB1']

for locus in Locus:
  if locus in ABC_loci:
    prefix = 'ABC'
  elif locus in DR_loci:
    prefix = 'DRB'
  else:
    prefix = 'DQ'

  # _quantity
  eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_quantity'] = pd.to_numeric(eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_quantity'], errors='coerce')
  eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_quantity'] = pd.to_numeric(eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_quantity'], errors='coerce')

  #eplet_dataframe[f'{prefix}_explanation_{locus}_quantity'] = eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_quantity'].fillna(0) + eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_quantity'].fillna(0)
  eplet_dataframe[f'{prefix}_explanation_{locus}_quantity'] = eplet_dataframe[[f'{prefix}_explanation_{locus}_allele_1_quantity', f'{prefix}_explanation_{locus}_allele_2_quantity']].sum(axis=1, min_count=1)

  # _details
  eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_details'].fillna(value='[]', inplace=True)
  eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_details'].fillna(value='[]', inplace=True)

  eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_details'] = eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_details'].apply(ast.literal_eval)
  eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_details'] = eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_details'].apply(ast.literal_eval)

  eplet_dataframe[f'{prefix}_explanation_{locus}_details'] = eplet_dataframe.apply(lambda row: row[f'{prefix}_explanation_{locus}_allele_1_details'] + row[f'{prefix}_explanation_{locus}_allele_2_details'], axis=1)

  # _explanation
  eplet_dataframe[f'{prefix}_explanation_{locus}_explanation'] = eplet_dataframe[f'{prefix}_explanation_{locus}_allele_1_explanation'].combine_first(eplet_dataframe[f'{prefix}_explanation_{locus}_allele_2_explanation'])

  # Drop allele columns
  eplet_dataframe.drop(columns=[
      f'{prefix}_explanation_{locus}_allele_1_quantity', f'{prefix}_explanation_{locus}_allele_1_details', f'{prefix}_explanation_{locus}_allele_1_explanation',
      f'{prefix}_explanation_{locus}_allele_2_quantity', f'{prefix}_explanation_{locus}_allele_2_details', f'{prefix}_explanation_{locus}_allele_2_explanation'
      ], inplace=True)

# Define a grouping order for the columns
group_order = ['PX', 'ALL', 'ABC', 'DRB', 'DQ', 'DP', 'MICA']

# Group and reorder columns
ordered_columns = []
for group in group_order:
    grouped_cols = [col for col in eplet_dataframe.columns if col.startswith(group)]
    ordered_columns.extend(grouped_cols)

# Reindex the DataFrame columns according to the defined order
eplet_dataframe = eplet_dataframe[ordered_columns]

# Processed by locus
eplet_dataframe.to_csv(output_file_path_PROC_loc, index=False, sep=',')

# Removing the prefix '{prefix}_explanation_' from column names
eplet_dataframe.columns = eplet_dataframe.columns.str.replace('ABC_explanation_', '').str.replace('DRB_explanation_', '').str.replace('DQ_explanation_', '')

# Cleaned version
eplet_dataframe.to_csv(output_file_path_PROC_loc_labels, index=False, sep=',')
