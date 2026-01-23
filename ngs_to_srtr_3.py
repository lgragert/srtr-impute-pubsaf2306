import os
import pandas as pd
import re
import pyard

baylor_file = 'Center_data_v2/Baylor_2017_10_2023_typing_for_Loren_update.tsv.gz'
# duke_file = 'Center_data_v2/Duke_2017_6_2023_typing_for_Loren_update.tsv.gz'
hopkins_file = 'Center_data_v2/Hopkins_main_table_2017_2024_for_Loren.tsv.gz'
penn_file = 'Center_data_v2/Penn_2017_2024_typing_for_Loren_update.tsv.gz'
readme_file = 'Center_data_v2/README_update.txt'
# stanford_file = 'Center_data_v2/Stanford_2017_6_2024_typing_for_Loren.tsv.gz'
ucsf_file = 'Center_data_v2/UCSF_2017_2024_typing_for_Loren.tsv.gz'
unc_file = 'Center_data_v2/UNC_2025_05_typing_for_Loren_update.tsv.gz'
uwisc_file = 'Center_data_v2/WIUW_2017_2023_typing_for_Loren_update.tsv.gz'

baylor_data = pd.read_csv(baylor_file, sep='\t', compression='gzip')
# duke_data = pd.read_csv(duke_file, sep='\t', compression='gzip')
hopkins_data = pd.read_csv(hopkins_file, sep='\t', compression='gzip')
penn_data = pd.read_csv(penn_file, sep='\t', compression='gzip')
# stanford_data = pd.read_csv(stanford_file, sep='\t', compression='gzip')
ucsf_data = pd.read_csv(ucsf_file, sep='\t', compression='gzip', encoding='ISO-8859-1') # or encoding='latin1'
unc_data = pd.read_csv(unc_file, sep='\t', compression='gzip')
uwisc_data = pd.read_csv(uwisc_file, sep='\t', compression='gzip')

print(f"Baylor data: {baylor_data.shape[0]} rows, {baylor_data.shape[1]} columns")
# print(f"Duke data: {duke_data.shape[0]} rows, {duke_data.shape[1]} columns")
print(f"Hopkins data: {hopkins_data.shape[0]} rows, {hopkins_data.shape[1]} columns")
print(f"Penn data: {penn_data.shape[0]} rows, {penn_data.shape[1]} columns")
# print(f"Stanford data: {stanford_data.shape[0]} rows, {stanford_data.shape[1]} columns")
print(f"UCSF data: {ucsf_data.shape[0]} rows, {ucsf_data.shape[1]} columns")
print(f"UNC data: {unc_data.shape[0]} rows, {unc_data.shape[1]} columns")
print(f"Wisconsin data: {uwisc_data.shape[0]} rows, {uwisc_data.shape[1]} columns")

print(f"Baylor columns: {baylor_data.columns.tolist()}")
# print(f"Duke columns: {duke_data.columns.tolist()}")
print(f"Hopkins columns: {hopkins_data.columns.tolist()}")
print(f"Penn columns: {penn_data.columns.tolist()}")
# print(f"Stanford columns: {stanford_data.columns.tolist()}")
print(f"UCSF columns: {ucsf_data.columns.tolist()}")
print(f"UNC columns: {unc_data.columns.tolist()}")
print(f"Wisconsin columns: {uwisc_data.columns.tolist()}")

baylor_NGS = baylor_data[baylor_data['typing_comment'].str.startswith('Generic high resolution', na=False)]
hopkins_NGS = hopkins_data[hopkins_data['typing_comment'].str.startswith('Generic high resolution', na=False)]
penn_NGS = penn_data[penn_data['typing_comment'].str.startswith('Generic high resolution', na=False)]
ucsf_NGS = ucsf_data[ucsf_data['typing_comment'].str.startswith('Generic high resolution', na=False)]
unc_NGS = unc_data[unc_data['GL_string'].notna()]
uwisc_NGS = uwisc_data[uwisc_data['typing_comment'].str.startswith('Generic high resolution', na=False)]

print(f"Baylor generic high resolution GL_string unique value: {baylor_NGS['GL_string'].dropna().unique()[0]}")
print(f"Hopkins generic high resolution GL_string unique value: {hopkins_NGS['GL_string'].dropna().unique()[0]}")
print(f"Penn generic high resolution GL_string unique value: {penn_NGS['GL_string'].dropna().unique()[0]}")
print(f"UCSF generic high resolution GL_string unique value: {ucsf_NGS['GL_string'].dropna().unique()[0]}")
print(f"UNC generic high resolution GL_string unique value: {unc_NGS['GL_string'].dropna().unique()[0]}")
print(f"Wisconsin generic high resolution GL_string unique value: {uwisc_NGS['GL_string'].dropna().unique()[0]}")

baylor_NGS_D = baylor_NGS[baylor_NGS['typing_ID'].str.startswith('d-')]
baylor_NGS_R = baylor_NGS[baylor_NGS['typing_ID'].str.startswith('r-')]

hopkins_NGS_D = hopkins_NGS[hopkins_NGS['typing_ID'].str.startswith('d-')]
hopkins_NGS_R = hopkins_NGS[hopkins_NGS['typing_ID'].str.startswith('r-')]

penn_NGS_D = penn_NGS[penn_NGS['typing_ID'].str.startswith('d-')]
penn_NGS_R = penn_NGS[penn_NGS['typing_ID'].str.startswith('r-')]

ucsf_NGS_D = ucsf_NGS[ucsf_NGS['typing_ID'].str.startswith('d-')]
ucsf_NGS_R = ucsf_NGS[ucsf_NGS['typing_ID'].str.startswith('r-')]

unc_NGS_D = unc_NGS[unc_NGS['typing_ID'].str.startswith('d-')]
unc_NGS_R = unc_NGS[unc_NGS['typing_ID'].str.startswith('r-')]

uwisc_NGS_D = uwisc_NGS[uwisc_NGS['typing_ID'].str.startswith('d-')]
uwisc_NGS_R = uwisc_NGS[uwisc_NGS['typing_ID'].str.startswith('r-')]

print(f"Baylor NGS Donor data: {baylor_NGS_D.shape[0]} rows, {baylor_NGS_D.shape[1]} columns")
print(f"Baylor NGS Recipient data: {baylor_NGS_R.shape[0]} rows, {baylor_NGS_R.shape[1]} columns")

print(f"Hopkins NGS Donor data: {hopkins_NGS_D.shape[0]} rows, {hopkins_NGS_D.shape[1]} columns")
print(f"Hopkins NGS Recipient data: {hopkins_NGS_R.shape[0]} rows, {hopkins_NGS_R.shape[1]} columns")

print(f"Penn NGS Donor data: {penn_NGS_D.shape[0]} rows, {penn_NGS_D.shape[1]} columns")
print(f"Penn NGS Recipient data: {penn_NGS_R.shape[0]} rows, {penn_NGS_R.shape[1]} columns")

print(f"UCSF NGS Donor data: {ucsf_NGS_D.shape[0]} rows, {ucsf_NGS_D.shape[1]} columns")
print(f"UCSF NGS Recipient data: {ucsf_NGS_R.shape[0]} rows, {ucsf_NGS_R.shape[1]} columns")

print(f"UNC NGS Donor data: {unc_NGS_D.shape[0]} rows, {unc_NGS_D.shape[1]} columns")
print(f"UNC NGS Recipient data: {unc_NGS_R.shape[0]} rows, {unc_NGS_R.shape[1]} columns")

print(f"Wisconsin NGS Donor data: {uwisc_NGS_D.shape[0]} rows, {uwisc_NGS_D.shape[1]} columns")
print(f"Wisconsin NGS Recipient data: {uwisc_NGS_R.shape[0]} rows, {uwisc_NGS_R.shape[1]} columns")

srtr_matrix = pd.read_csv(srtr_file, sep='\t', compression='gzip')
print(f'SRTR matrix data: {srtr_matrix.shape[0]} rows, {srtr_matrix.shape[1]} columns')
print(f'SRTR matrix columns: {srtr_matrix.columns.tolist()}')

def parse_hla_typing(df, drtype):
    """
    Parse HLA GL_string format and create individual allele columns
    
    Parameters:
    df: DataFrame containing GL_string column
    drtype: Either 'DONOR' or 'RECIP' to prefix column names
    
    Returns:
    DataFrame with new allele columns added
    """
    df = df.copy()
    
    # Initialize all locus columns
    loci_cols = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    for locus in loci_cols:
        df[f'{drtype}_{locus}'] = ''
        df[f'{drtype}_{locus}_1'] = ''
        df[f'{drtype}_{locus}_2'] = ''
    
    for idx, row in df.iterrows():

        # if pd.isna(row['GL_string']) or row['GL_string'] == '':
        #     continue
            
        # Split GL_string by ^ to get individual locus pieces
        pieces = row['GL_string'].split('^')
        
        for piece in pieces:
            piece = piece.strip()
            if not piece:
                continue
                
            locus_assigned = False
            if piece.startswith('HLA-A*'): # TODO - UNC-donor-df has some typings without the *-prefix (doesn't get assigned to any locus)
                locus = 'A'
                locus_assigned = True
            elif piece.startswith('HLA-B*'):
                locus = 'B'
                locus_assigned = True
            elif piece.startswith('HLA-C*'):
                locus = 'C'
                locus_assigned = True
            elif piece.startswith('HLA-DRB1*'):
                locus = 'DRB1'
                locus_assigned = True
            elif piece.startswith('HLA-DRB3*') or piece.startswith('HLA-DRB4*') or piece.startswith('HLA-DRB5*'):
                locus = 'DRB345'
                locus_assigned = True
            elif piece.startswith('HLA-DQA1*'):
                locus = 'DQA1'
                locus_assigned = True
            elif piece.startswith('HLA-DQB1*'):
                locus = 'DQB1'
                locus_assigned = True
            elif piece.startswith('HLA-DPA1*'):
                locus = 'DPA1'
                locus_assigned = True
            elif piece.startswith('HLA-DPB1*'):
                locus = 'DPB1'
                locus_assigned = True
            
            if locus_assigned:
                # Assign to intermediate column
                df.at[idx, f'{drtype}_{locus}'] = piece

                alleles = piece.split('+')
                alleles = [allele.strip().replace('HLA-', '') for allele in alleles]
                
                if len(alleles) >= 1 and alleles[0]:
                    df.at[idx, f'{drtype}_{locus}_1'] = alleles[0]
                if len(alleles) >= 2 and alleles[1]:
                    df.at[idx, f'{drtype}_{locus}_2'] = alleles[1]
                elif len(alleles) == 1:  # Homozygous case
                    if locus == 'DRB345':
                        # Special handling for DRB345 - set DRBX*NNNN for empty DRB345 when DRB1 is *01, *08, or *10
                        drb1_allele_1 = df.at[idx, f'{drtype}_DRB1_1']
                        drb1_allele_2 = df.at[idx, f'{drtype}_DRB1_2']
                        if re.match(r'DRB1\*(01|08|10)', drb1_allele_1) or re.match(r'DRB1\*(01|08|10)', drb1_allele_2):
                            df.at[idx, f'{drtype}_{locus}_2'] = 'DRBX*NNNN'
                        else:
                            df.at[idx, f'{drtype}_{locus}_2'] = alleles[0]
                    else:
                        df.at[idx, f'{drtype}_{locus}_2'] = alleles[0]

        # If any DRB345 locus column is completely empty, set both locus alleles to 'DRBX*NNNN'
        if df.at[idx, f'{drtype}_DRB345'] == '':
            df.at[idx, f'{drtype}_DRB345_1'] = 'DRBX*NNNN'
            df.at[idx, f'{drtype}_DRB345_2'] = 'DRBX*NNNN'
    
    return df


hla_baylor_d = parse_hla_typing(baylor_NGS_D, 'DONOR')
hla_baylor_r = parse_hla_typing(baylor_NGS_R, 'RECIP')

hla_hopkins_d = parse_hla_typing(hopkins_NGS_D, 'DONOR')
hla_hopkins_r = parse_hla_typing(hopkins_NGS_R, 'RECIP')

hla_penn_d = parse_hla_typing(penn_NGS_D, 'DONOR')
hla_penn_r = parse_hla_typing(penn_NGS_R, 'RECIP')

hla_ucsf_d = parse_hla_typing(ucsf_NGS_D, 'DONOR')
hla_ucsf_r = parse_hla_typing(ucsf_NGS_R, 'RECIP')

hla_unc_d = parse_hla_typing(unc_NGS_D, 'DONOR')
hla_unc_r = parse_hla_typing(unc_NGS_R, 'RECIP')

hla_uwisc_d = parse_hla_typing(uwisc_NGS_D, 'DONOR')
hla_uwisc_r = parse_hla_typing(uwisc_NGS_R, 'RECIP')

if not hla_baylor_d.empty:
    print("Baylor Donor - First row GL_string:")
    print(baylor_NGS_D['GL_string'].iloc[0])
    print("Parsed alleles:")
    allele_cols = [col for col in hla_baylor_d.columns if col.startswith('DONOR_') and ('_1' in col or '_2' in col)]
    for col in sorted(allele_cols):
        print(f"  {col}: {hla_baylor_d[col].iloc[0]}")

print(f"\nColumn check - Baylor Donor has {len([col for col in hla_baylor_d.columns if col.startswith('DONOR_')])} DONOR columns")

for center_name, df_donor, df_recip in [
    ('Baylor', hla_baylor_d, hla_baylor_r),
    ('Hopkins', hla_hopkins_d, hla_hopkins_r),
    ('Penn', hla_penn_d, hla_penn_r),
    ('UCSF', hla_ucsf_d, hla_ucsf_r),
    ('UNC', hla_unc_d, hla_unc_r),
    ('Wisconsin', hla_uwisc_d, hla_uwisc_r)
]:
    if not df_donor.empty:
        donor_cols = len([col for col in df_donor.columns if col.startswith('DONOR_')])
        print(f"{center_name} Donor: {donor_cols} DONOR columns")
    if not df_recip.empty:
        recip_cols = len([col for col in df_recip.columns if col.startswith('RECIP_')])
        print(f"{center_name} Recipient: {recip_cols} RECIP columns")

def clean_hla_data(df, drtype):
    """
    Clean HLA typing data by removing rows with incomplete or invalid alleles
    
    Parameters:
    df: DataFrame containing parsed HLA allele columns
    drtype: Either 'DONOR' or 'RECIP' to identify column prefix
    
    Returns:
    DataFrame with invalid alleles removed
    """
    df_clean = df.copy()
    loci_cols = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    
    print(f"Before cleaning, total rows: {df_clean.shape[0]}")

    # Drop rows with empty allele fields
    for locus in loci_cols:
        col1 = f'{drtype}_{locus}_1'
        col2 = f'{drtype}_{locus}_2'
        mask = df_clean[col1].isna() | df_clean[col1].eq('') | df_clean[col2].isna() | df_clean[col2].eq('')
        print(f"List of rows with empty alleles in {col1} or {col2}:")
        print(df_clean[mask][['typing_ID', col1, col2]])
        df_clean = df_clean[~mask]
        print(f"After dropping empty alleles in {col1} or {col2}, remaining rows: {df_clean.shape[0]}")

    # Drop rows with alleles with invalid alleles
    for locus in loci_cols:
        if locus == 'DRB345':
            allele_prefixes = ['DRB3', 'DRB4', 'DRB5']
        else:
            allele_prefixes = [locus]

        for allele_prefix in allele_prefixes:
            col1 = f'{drtype}_{locus}_1'
            col2 = f'{drtype}_{locus}_2'

            # Drop rows with alleles with values NOT containing any colons in it (e.g., A*02, A*0301) - EXCEPT for DRBX*NNNN in DRB345 columns
            if locus == 'DRB345':
                # For DRB345, keep DRBX*NNNN values but drop others without colons
                mask = (~df_clean[col1].str.contains(rf':', na=False) & (df_clean[col1] != 'DRBX*NNNN')) | (~df_clean[col2].str.contains(rf':', na=False) & (df_clean[col2] != 'DRBX*NNNN'))
            else:
                # For other loci, drop any values without colons
                mask = ~df_clean[col1].str.contains(rf':', na=False) | ~df_clean[col2].str.contains(rf':', na=False)
            print(f"List of alleles without colons in {col1} or {col2}:")
            print(df_clean[mask][['typing_ID', col1, col2]])
            df_clean = df_clean[~mask]
            print(f"After dropping alleles without colons in {col1} or {col2}, remaining rows: {df_clean.shape[0]}")

            # Drop rows NOT containing any asterisks in it (e.g., A11)
            mask = ~df_clean[col1].str.contains(rf'\*', na=False) | ~df_clean[col2].str.contains(rf'\*', na=False)
            print(f"List of alleles without asterisks in {col1} or {col2}:")
            print(df_clean[mask][['typing_ID', col1, col2]])
            df_clean = df_clean[~mask]
            print(f"After dropping alleles without asterisks in {col1} or {col2}, remaining rows: {df_clean.shape[0]}")

            # Drop rows with single-field alleles containing 2 or 3 digits or letters (e.g., A*31 or DPA1*XX)
            mask = df_clean[col1].str.match(rf'{allele_prefix}\*[A-Za-z0-9]{{2,3}}$', na=False) | df_clean[col2].str.match(rf'{allele_prefix}\*[A-Za-z0-9]{{2,3}}$', na=False)
            print(f"List of single-field alleles in {col1} or {col2}:")
            print(df_clean[mask][['typing_ID', col1, col2]])
            df_clean = df_clean[~mask]
            print(f"After dropping single-field alleles in {col1} or {col2}, remaining rows: {df_clean.shape[0]}")

            # Drop rows with two- or three-digit incomplete alleles (e.g., A*31: or B*55: or DPB1*463:)
            mask = df_clean[col1].str.match(rf'{allele_prefix}\*[A-Za-z0-9]{{2,3}}:$', na=False) | df_clean[col2].str.match(rf'{allele_prefix}\*[A-Za-z0-9]{{2,3}}:$', na=False)
            print(f"List of incomplete alleles in {col1} or {col2}:")
            print(df_clean[mask][['typing_ID', col1, col2]])
            df_clean = df_clean[~mask]
            print(f"After dropping incomplete alleles in {col1} or {col2}, remaining rows: {df_clean.shape[0]}")


            # Truncate trailing colons from two-field alleles (e.g., DPA1*01:03:)
            # If the column contains "{allele_prefix}*{two-digit}:{two-digit}:", print those rows
            mask = df_clean[col1].str.match(rf'{allele_prefix}\*\d{{2}}:\d{{2}}:$', na=False) | df_clean[col2].str.match(rf'{allele_prefix}\*\d{{2}}:\d{{2}}:$', na=False)
            print(f"List of alleles with trailing colon in {col1} or {col2}:")
            print(df_clean[mask][['typing_ID', col1, col2]])
            df_clean[col1] = df_clean[col1].str.replace(rf'({allele_prefix}\*\d{{2}}:\d{{2}}):$', r'\1', regex=True)
            df_clean[col2] = df_clean[col2].str.replace(rf'({allele_prefix}\*\d{{2}}:\d{{2}}):$', r'\1', regex=True)
            print(f"After removing trailing colon in {col1} or {col2}, sample rows:")
            print(df_clean[mask][['typing_ID', col1, col2]])

    print(f"After cleaning, remaining rows: {df_clean.shape[0]}")

    return df_clean

all_center_dfs = [
    ('Baylor', hla_baylor_d, hla_baylor_r),
    ('Hopkins', hla_hopkins_d, hla_hopkins_r),
    ('Penn', hla_penn_d, hla_penn_r),
    ('UCSF', hla_ucsf_d, hla_ucsf_r),
    ('UNC', hla_unc_d, hla_unc_r),
    ('WISC', hla_uwisc_d, hla_uwisc_r)
]

cleaned_center_dfs = {}
for center_name, df_donor, df_recip in all_center_dfs:
    print(f"\n--- Processing {center_name} ---")
    
    if not df_donor.empty:
        print(f"Cleaning {center_name} Donor Data")
        cleaned_donor = clean_hla_data(df_donor, 'DONOR')
        cleaned_center_dfs[f'clean_{center_name.lower()}_donor'] = cleaned_donor
    
    if not df_recip.empty:
        print(f"Cleaning {center_name} Recipient Data")
        cleaned_recip = clean_hla_data(df_recip, 'RECIP')
        cleaned_center_dfs[f'clean_{center_name.lower()}_recip'] = cleaned_recip

clean_Baylor_donor = cleaned_center_dfs['clean_baylor_donor']
clean_Baylor_recip = cleaned_center_dfs['clean_baylor_recip']

clean_Hopkins_donor = cleaned_center_dfs['clean_hopkins_donor']
clean_Hopkins_recip = cleaned_center_dfs['clean_hopkins_recip']

clean_Penn_donor = cleaned_center_dfs['clean_penn_donor']
clean_Penn_recip = cleaned_center_dfs['clean_penn_recip']

clean_UCSF_donor = cleaned_center_dfs['clean_ucsf_donor']
clean_UCSF_recip = cleaned_center_dfs['clean_ucsf_recip']

clean_UNC_donor = cleaned_center_dfs['clean_unc_donor']
clean_UNC_recip = cleaned_center_dfs['clean_unc_recip']

clean_WISC_donor = cleaned_center_dfs['clean_wisc_donor']
clean_WISC_recip = cleaned_center_dfs['clean_wisc_recip']

ard = pyard.init()

def pyard_u2_redux(df, drtype):
    """
    Reduce alleles in the dataframe using pyard
    
    Parameters:
    df: DataFrame containing allele columns
    drtype: Either 'DONOR' or 'RECIP' to prefix column names
    
    Returns:
    DataFrame with reduced allele columns
    """
    df_u2 = df.copy()
    
    allele_cols = [col for col in df_u2.columns if col.startswith(f'{drtype}_') and ('_1' in col or '_2' in col)]
    print(f"Allele columns to reduce with pyard: {allele_cols}")

    for col in allele_cols:
        df_u2[col] = df_u2[col].apply(lambda x: ard.redux(x, 'U2') if x != 'DRBX*NNNN' else x)
    
    return df_u2


pyard_center_dfs = {}
for center_name, df_donor, df_recip in all_center_dfs:
    print(f"\n--- Applying pyard reduction for {center_name} ---")
    
    clean_donor = cleaned_center_dfs.get(f'clean_{center_name.lower()}_donor')
    clean_recip = cleaned_center_dfs.get(f'clean_{center_name.lower()}_recip')
    
    if clean_donor is not None:
        print(f"Reducing alleles for {center_name} Donor Data")
        reduced_donor = pyard_u2_redux(clean_donor, 'DONOR')
        pyard_center_dfs[f'pyard_{center_name.lower()}_donor'] = reduced_donor
    
    if clean_recip is not None:
        print(f"Reducing alleles for {center_name} Recipient Data")
        reduced_recip = pyard_u2_redux(clean_recip, 'RECIP')
        pyard_center_dfs[f'pyard_{center_name.lower()}_recip'] = reduced_recip

pyard_Baylor_donor = pyard_center_dfs['pyard_baylor_donor']
pyard_Baylor_recip = pyard_center_dfs['pyard_baylor_recip']

pyard_Hopkins_donor = pyard_center_dfs['pyard_hopkins_donor']
pyard_Hopkins_recip = pyard_center_dfs['pyard_hopkins_recip']

pyard_Penn_donor = pyard_center_dfs['pyard_penn_donor']
pyard_Penn_recip = pyard_center_dfs['pyard_penn_recip']

pyard_UCSF_donor = pyard_center_dfs['pyard_ucsf_donor']
pyard_UCSF_recip = pyard_center_dfs['pyard_ucsf_recip']

pyard_UNC_donor = pyard_center_dfs['pyard_unc_donor']
pyard_UNC_recip = pyard_center_dfs['pyard_unc_recip']

pyard_WISC_donor = pyard_center_dfs['pyard_wisc_donor']
pyard_WISC_recip = pyard_center_dfs['pyard_wisc_recip']

for key, df in pyard_center_dfs.items():

    output_file = f"{key}.csv"
    df.to_csv(output_file, index=False)
