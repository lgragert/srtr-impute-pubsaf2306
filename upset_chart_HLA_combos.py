
import pandas as pd
import itertools

# Create an input file format that has a list of HLA category combos
# Column headers as HLA categories
# Rows will have indicator variables (0,1)
# Comment field in the input to explain the rationale for looking at the combo
# DRB1_AgMM | DQA1_AgMM  # Example
#    0            0
#    1            0
#    0            1
#    1            1


# There are some combos that are not realistic in the numbers one: DRB1_0AgMM cannot exist at the same time as DRB1_1AgMM or 2AgMM, etc.
def is_valid_combo(row):

    loci_list = ['DRB1', 'DQA1', 'DQB1']  # Add to this list if you need to add more genes

    for locus in loci_list:
        # Cannot have these for Antigen MM
        if row[locus + '_0AgMM'] == 1 and row[locus + '_1AgMM'] == 1:
            return False
        if row[locus + '_0AgMM'] == 1 and row[locus + '_2AgMM'] == 1:
            return False
        if row[locus + '_1AgMM'] == 1 and row[locus + '_2AgMM'] == 1:
            return False
        if row[locus + '_0AgMM'] == 0 and row[locus + '_1AgMM'] == 0 and row[locus + '_2AgMM'] == 0:  # one of the categories has to exist
            return False

        # cannot have these for Allele MM
        if row[locus + '_0AlleleMM'] == 1 and row[locus + '_1AlleleMM'] == 1:
            return False
        if row[locus + '_0AlleleMM'] == 1 and row[locus + '_2AlleleMM'] == 1:
            return False
        if row[locus + '_1AlleleMM'] == 1 and row[locus + '_2AlleleMM'] == 1:
            return False
        if row[locus + '_0AlleleMM'] == 0 and row[locus + '_1AlleleMM'] == 0 and row[locus + '_2AlleleMM'] == 0:  # one of the categories has to exist
            return False

    return True


# We want to look at Class II (minus DP) AlleleMM and AgMM
HLA_headers = ['DRB1_AgMM', 'DRB1_AlleleMM', 'DQA1_AgMM', 'DQA1_AlleleMM', 'DQB1_AgMM', 'DQB1_AlleleMM']
HLA_header_nums = ['DRB1_0AgMM', 'DRB1_1AgMM', 'DRB1_2AgMM', 'DRB1_0AlleleMM', 'DRB1_1AlleleMM', 'DRB1_2AlleleMM',
                   'DQA1_0AgMM', 'DQA1_1AgMM', 'DQA1_2AgMM', 'DQA1_0AlleleMM', 'DQA1_1AlleleMM', 'DQA1_2AlleleMM',
                   'DQB1_0AgMM', 'DQB1_1AgMM', 'DQB1_2AgMM', 'DQB1_0AlleleMM', 'DQB1_1AlleleMM', 'DQB1_2AlleleMM']

# Generate all possible combinations of indicator variables (0, 1) for the given column
HLA_combinations = list(itertools.product([0,1], repeat=len(HLA_headers)))
HLA_combinations_nums = list(itertools.product([0,1], repeat=len(HLA_header_nums)))

# Create the DataFrame from these combinations
HLA_combos = pd.DataFrame(HLA_combinations, columns=HLA_headers)
HLA_combos_num = pd.DataFrame(HLA_combinations_nums, columns=HLA_header_nums)


# Apply the filter function to the DataFrame to only have valid combinations
valid_combo_nums = HLA_combos_num.apply(is_valid_combo, axis=1)
filtered_combos_nums = HLA_combos_num[valid_combo_nums]
filtered_combos_nums = filtered_combos_nums.reset_index(drop=True)

with pd.ExcelWriter('Upset_Chart_HLA_Combos.xlsx') as writer:
    HLA_combos.to_excel(writer, sheet_name='General_MM', index=False)
    filtered_combos_nums.to_excel(writer, sheet_name='Numbered_MM', index=False)
