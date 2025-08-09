
# Create an input file with columns we want to compare in the upset chart
# Column headers: MM_pos_# (0 or 1) vs MM_Ag_# directional present in donor not recip (0 or 1)
# Rows will have indicator variables (0,1)

import os
import pandas as pd
import pyard



ard = pyard.init()
print("Initialized pyard with the latest version of IPD-IMGT/HLA database.")


SRTR_matrix = pd.read_csv('SRTR_AA_MM_9loc_matrix_1.txt.gz', sep='\t', compression='gzip')
print(f"Read {len(SRTR_matrix)} rows from the SRTR matrix file.")

# %%

SRTR_matrix = SRTR_matrix[SRTR_matrix['ORG_TY'] == 'KI: Kidney']
print(f"Filtered to {len(SRTR_matrix)} rows with ORG_TY 'KI: Kidney'.")

SRTR_matrix = SRTR_matrix[SRTR_matrix['DON_TY'] == 'C']
print(f"Filtered to {len(SRTR_matrix)} rows with DON_TY 'C'.")

# Handle string values like '35-49' by converting to numeric or filtering them out
SRTR_matrix['REC_AGE_AT_TX'] = pd.to_numeric(SRTR_matrix['REC_AGE_AT_TX'], errors='coerce')
SRTR_matrix = SRTR_matrix[SRTR_matrix['REC_AGE_AT_TX'] >= 18]
print(f"Filtered to {len(SRTR_matrix)} rows with REC_AGE_AT_TX >= 18.")

SRTR_matrix = SRTR_matrix[SRTR_matrix['DON_AGE'] >= 9.0]
print(f"Filtered to {len(SRTR_matrix)} rows with DON_AGE >= 9.0.")

# %%

# AA positions of interest (columns)
# Column datatype is int64
AA_POSITIONS = ['MM_A_17', 'MM_A_73', 'MM_A_76', 'MM_B_103', 'MM_C_114']


for pos in AA_POSITIONS:
    col_name = f'{pos}'
    new_col_name = f'{pos}_binary'

    # Check if any rows with NaN in the original column
    print(SRTR_matrix[[f'{pos}']].isna().sum())

    # Create a new column with binary values based on the original column
    SRTR_matrix[new_col_name] = SRTR_matrix[col_name].apply(lambda x: 1 if x in [1, 2] else 0)
    print(f"Created column {new_col_name} with values 0 or 1 based on {pos} values.")

    # Check if any NaN values in the new column
    print(SRTR_matrix[[new_col_name]].isna().sum())
print(f"Created new columns for AA positions: {AA_POSITIONS}")

# %%

# Two-field columns for donor and recipient
Ag_COLUMNS = ['DONOR_A_1', 'DONOR_A_2', 'RECIP_A_1', 'RECIP_A_2', 
              'DONOR_B_1', 'DONOR_B_2', 'RECIP_B_1', 'RECIP_B_2', 
              'DONOR_C_1', 'DONOR_C_2', 'RECIP_C_1', 'RECIP_C_2']

# Create new columns for First-field Antigen
# Take the original column value and reduce it to first-field serological level to assign to new column

for col in Ag_COLUMNS:
    print(SRTR_matrix[[col]].isna().sum())

    ff_col_name = f'{col}_first'
    SRTR_matrix[ff_col_name] = SRTR_matrix[col].apply(lambda x: ard.redux(x, 'S') if pd.notna(x) else None)

    print(SRTR_matrix[[ff_col_name]].isna().sum())
    print(f"Created column {ff_col_name} with first-field serological level values based on {col} values.")

# %%

# check column values
for col in Ag_COLUMNS:
    ff_col_name = f'{col}_first'
    unique_values = SRTR_matrix[ff_col_name].unique()
    print(f"Unique values in {ff_col_name}: {unique_values}")

# %%

# Looking at antigens present in donor but not in recipient - from KH bubble plots for ABC loci
Ag_COLUMNS_binary = ['A3,A33',
                     'A30,A2',
                     'A30,A29',
                     'A30,A1',
                     'A33,A68',
                     'A34,A24',
                     'A74,A3',
                     'A74,A24',
                     'A74,A1',
                     'A74,A68',
                     'B44,B8',
                     'B45,B8',
                     'B53,B8',
                     'B53,B44',
                     'B55,B27',
                     'B57,B41',
                     'B62,B49',
                     'B65,B49',
                     'B72,B8',
                     'B65,B8',
                     'Cw4,Cw7',
                     'Cw6,Cw7',
                     'Cw8,Cw7',
                     'C12,Cw3',
                     'C16,Cw7',
                     'C17,Cw6',
                     'C18,Cw7',
                     'C18,Cw8',
                     'C18,Cw4',
                     'C18,Cw6']

# Initialize new columns in the SRTR matrix for binary values
# Each column will indicate if the donor has the antigen and the recipient does not
for col in Ag_COLUMNS_binary:
    print(f"Processing column {col}...")

    # Initialize the binary column with 0
    SRTR_matrix[col] = 0
    print(f"Initialized column {col} with 0.")

    # Indicate which locus variable by the first letter of the column name
    loc = col[0]  # 'A', 'B', or 'C'
    print(f"Processing column {col} for locus {loc}...")

    # Indicate which antigen is being checked for presence in donor and recipient
    # e.g. for 'A34,A24', donor has A34 and recipient has A24
    don_ag, recip_ag = col.split(',')
    print(f"Donor antigen: {don_ag}, Recipient antigen: {recip_ag}")


    donor_has_ag = ((SRTR_matrix[f'DONOR_{loc}_1_first'].str.contains(don_ag, na=False)) |
                    (SRTR_matrix[f'DONOR_{loc}_2_first'].str.contains(don_ag, na=False)))
    recip_has_ag = ((SRTR_matrix[f'RECIP_{loc}_1_first'].str.contains(recip_ag, na=False)) |
                    (SRTR_matrix[f'RECIP_{loc}_2_first'].str.contains(recip_ag, na=False)))
    recip_no_ag = ((~SRTR_matrix[f'RECIP_{loc}_1_first'].str.contains(don_ag, na=False)) &
                   (~SRTR_matrix[f'RECIP_{loc}_2_first'].str.contains(don_ag, na=False)))
    SRTR_matrix[col] = (donor_has_ag & recip_has_ag & recip_no_ag).astype(int)
    

    print(f"Created column {col} with binary values based on first-field values.")


# %%

upset_columns = ['TX_ID', 'PX_ID'] + \
                [f'{pos}_binary' for pos in AA_POSITIONS] + \
                [f'{col}_first' for col in Ag_COLUMNS] + \
                Ag_COLUMNS_binary
print(f"Selected columns for the upset chart: {upset_columns}")

SRTR_upset_data = SRTR_matrix[upset_columns]
print(f"Subsetted DataFrame to {len(SRTR_upset_data)} rows and {len(SRTR_upset_data.columns)} columns for the upset chart.")
print(SRTR_upset_data.head())

output_filename = "SRTR_AA_MM_9loc_matrix_1_upset_all.csv"
SRTR_upset_data.to_csv(output_filename, index=False)
print(f"Saved the updated SRTR matrix with upset categories to {output_filename}.")

# %%


SRTR_UPSET = pd.read_csv('SRTR_AA_MM_9loc_matrix_1_upset.csv')
print(f"Read {len(SRTR_UPSET)} rows from the SRTR matrix upset file.")

# %%

# Create a list of AgMM combinations from KH bubble plots
agmm_combinations = ['A3,A33','A30,A2','A30,A29','A30,A1','A33,A68','A34,A24','A74,A3','A74,A24','A74,A1','A74,A68',
                     'B44,B8','B45,B8','B53,B8','B53,B44','B55,B27','B57,B41','B62,B49','B65,B49','B72,B8','B65,B8',
                     'Cw4,Cw7','Cw6,Cw7','Cw8,Cw7','C12,Cw3','C16,Cw7','C17,Cw6','C18,Cw7','C18,Cw8','C18,Cw4','C18,Cw6']

# Initialize a dictionary to store the counts and fractions
agmm_counts = {}
for combo in agmm_combinations:
    locus = combo[0]  # 'A', 'B', or 'C'
    donor_ag, recip_ag = combo.split(',')
    col_name = f'{donor_ag}_{recip_ag}'
    SRTR_UPSET[col_name] = 0
    print(f"Initialized column {col_name} with 0.")


    donor_has_donor_ag = ((SRTR_UPSET[f'DONOR_{locus}_1_first'].str.contains(donor_ag, na=False)) | 
                          (SRTR_UPSET[f'DONOR_{locus}_2_first'].str.contains(donor_ag, na=False)))
    recip_has_recip_ag = ((SRTR_UPSET[f'RECIP_{locus}_1_first'].str.contains(recip_ag, na=False)) | 
                          (SRTR_UPSET[f'RECIP_{locus}_2_first'].str.contains(recip_ag, na=False)))
    recip_no_donor_ag = ((~SRTR_UPSET[f'RECIP_{locus}_1_first'].str.contains(donor_ag, na=False)) & 
                         (~SRTR_UPSET[f'RECIP_{locus}_2_first'].str.contains(donor_ag, na=False)))
    SRTR_UPSET[col_name] = (donor_has_donor_ag & recip_has_recip_ag & recip_no_donor_ag).astype(int)

    print(f"Created column {col_name} with binary values based on first-field values for {combo}.")

    if col_name in SRTR_UPSET.columns:

        # Count the number of instances where the combination is present
        count = SRTR_UPSET[col_name].sum()
        print(f"Count for {col_name}: {count}")

        # Calculate the fraction of total instances
        fraction = count / len(SRTR_UPSET) if len(SRTR_UPSET) > 0 else 0
        print(f"Fraction for {col_name}: {fraction:.4f}")

        # Store the count and fraction in the dictionary
        agmm_counts[combo] = {'count': count, 'fraction': fraction}

    else:
        print(f"Column {col_name} not found in SRTR_UPSET DataFrame.")

for combo, data in agmm_counts.items():
    print(f"Combination {combo}: Count = {data['count']}, Fraction = {data['fraction']:.4f}")

agmm_counts_df = pd.DataFrame.from_dict(agmm_counts, orient='index')
agmm_counts_df.reset_index(inplace=True)
agmm_counts_df.columns = ['AgMM_Combo', 'Count', 'Fraction']
output_agmm_filename = "SRTR_AgMM_combinations_table.csv"
agmm_counts_df.to_csv(output_agmm_filename, index=False)
print(f"Saved the AgMM combinations counts and fractions to {output_agmm_filename}.")

