import pandas as pd
import hlagenie

# Create a file that has all the AA positions and AAMMs across all pairs
genie = hlagenie.init('3510')

# Read in SRTR matrix files 1-10
for num in range(1,2):
    srtr_filename = 'SRTR_AA_MM_9loc_matrix_' + str(num) + '.txt.gz'
    srtr_matrix = pd.read_csv(srtr_filename, sep='\t', header=0, compression='gzip')
    print(f'Processing file: {srtr_filename}')

    # Only need columns where it shows the AAMM and the haplotypes
    mm_cols = srtr_matrix.columns[srtr_matrix.columns.str.startswith('MM_')].tolist()
    patient_cols = ['TX_ID', 'PX_ID', 'HAPPAIR_RECIP', 'HAPPAIR_DONOR']
    select_cols = patient_cols + mm_cols
    srtr_matrix1 = srtr_matrix[select_cols]

    # Break up happair columns into hap1 and hap2 into the individual alleles
    srtr_matrix1[['HAPPAIR_RECIP_1', 'HAPPAIR_RECIP_2']] = srtr_matrix1['HAPPAIR_RECIP'].str.split('+', expand=True)
    srtr_matrix1[['A_RECIP_1', 'C_RECIP_1', 'B_RECIP_1', 'DRB345_RECIP_1', 'DRB1_RECIP_1', 'DQA1_RECIP_1', 'DQB1_RECIP_1',
              'DPA1_RECIP_1', 'DPB1_RECIP_1']] = srtr_matrix1['HAPPAIR_RECIP_1'].str.split('~', expand=True)
    srtr_matrix1[['A_RECIP_2', 'C_RECIP_2', 'B_RECIP_2', 'DRB345_RECIP_2', 'DRB1_RECIP_2', 'DQA1_RECIP_2', 'DQB1_RECIP_2',
              'DPA1_RECIP_2', 'DPB1_RECIP_2']] = srtr_matrix1['HAPPAIR_RECIP_2'].str.split('~', expand=True)

    srtr_matrix1[['HAPPAIR_DONOR_1', 'HAPPAIR_DONOR_2']] = srtr_matrix1['HAPPAIR_DONOR'].str.split('+', expand=True)
    srtr_matrix1[['A_DONOR_1', 'C_DONOR_1', 'B_DONOR_1', 'DRB345_DONOR_1', 'DRB1_DONOR_1', 'DQA1_DONOR_1', 'DQB1_DONOR_1',
              'DPA1_DONOR_1', 'DPB1_DONOR_1']] = srtr_matrix1['HAPPAIR_DONOR_1'].str.split('~', expand=True)
    srtr_matrix1[['A_DONOR_2', 'C_DONOR_2', 'B_DONOR_2', 'DRB345_DONOR_2', 'DRB1_DONOR_2', 'DQA1_DONOR_2', 'DQB1_DONOR_2',
              'DPA1_DONOR_2', 'DPB1_DONOR_2']] = srtr_matrix1['HAPPAIR_DONOR_2'].str.split('~', expand=True)

    # Get AAs from each position we have a MM
    hla_loc = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
    for hla in hla_loc:
        # Get the MM columns for this HLA
        mm_cols = srtr_matrix1.columns[srtr_matrix1.columns.str.startswith('MM_' + hla) & ~srtr_matrix1.columns.str.endswith('COUNT')].tolist()
        mm_pos = len(mm_cols)

        # Get the AAs from the MM columns
        print(f'Getting AAs for {hla} in file {num}')
        for pos in range(1, mm_pos + 1):
            aa_col = 'AA_' + hla + '_' + str(pos)

            # Applies genie to get the AAs
            srtr_matrix1[aa_col + '_RECIP_1'] = srtr_matrix1[hla + '_RECIP_1'].apply(lambda x: genie.getAA(x, pos) if x != 'DRBX*NNNN' else '-')
            srtr_matrix1[aa_col + '_RECIP_2'] = srtr_matrix1[hla + '_RECIP_2'].apply(lambda x: genie.getAA(x, pos) if x != 'DRBX*NNNN' else '-')
            srtr_matrix1[aa_col + '_DONOR_1'] = srtr_matrix1[hla + '_DONOR_1'].apply(lambda x: genie.getAA(x, pos) if x != 'DRBX*NNNN' else '-')
            srtr_matrix1[aa_col + '_DONOR_2'] = srtr_matrix1[hla + '_DONOR_2'].apply(lambda x: genie.getAA(x, pos) if x != 'DRBX*NNNN' else '-')

        # Clean up the DataFrame to only have the relevant columns
        hla_cols = ['TX_ID', 'PX_ID', 'HAPPAIR_RECIP', 'HAPPAIR_DONOR']
        aa_cols = srtr_matrix1.columns[srtr_matrix1.columns.str.startswith('AA_' + hla)].tolist()
        mm_cols = srtr_matrix1.columns[srtr_matrix1.columns.str.startswith('MM_' + hla) & ~srtr_matrix1.columns.str.endswith('COUNT')].tolist()
        select_cols = hla_cols + aa_cols + mm_cols
        srtr_matrix2 = srtr_matrix1[select_cols]

        # Put columns in the order that Malek wants, AA_RECIP1/2, AA_DONOR1/2, MM cols
        ordered_cols = []
        mm_cols = srtr_matrix1.columns[
            srtr_matrix1.columns.str.startswith('MM_' + hla) & ~srtr_matrix1.columns.str.endswith('COUNT')].tolist()
        mm_pos = len(mm_cols)
        for pos in range(1, mm_pos + 1):
            for suffix in ['RECIP_1', 'RECIP_2', 'DONOR_1', 'DONOR_2']:
                col_name = f'AA_{hla}_{pos}_{suffix}'
                if col_name in srtr_matrix2.columns:
                    ordered_cols.append(col_name)
            mm_col = f'MM_{hla}_{pos}'
            if mm_col in srtr_matrix2.columns:
                ordered_cols.append(mm_col)

        # Add patient columns at the start
        patient_cols = ['TX_ID', 'PX_ID', 'HAPPAIR_RECIP', 'HAPPAIR_DONOR']
        final_cols = patient_cols + ordered_cols

        # Reindex DataFrame
        srtr_matrix3 = srtr_matrix2[final_cols]
        print(f'Writing file: SRTR_AA_MM_9loc_with_HLA_{hla}_AA_{str(num)}.csv.gz')

        srtr_matrix3.to_csv(f'SRTR_matrix_{str(num)}_AA_HLA_{hla}.txt.gz', sep='\t', index=False, header=True, compression='gzip')
        print(f'Wrote file: SRTR_matrix_{str(num)}_AA_HLA_{hla}.txt.gz')
