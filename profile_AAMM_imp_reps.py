
import pandas as pd

# Profile the AAMM frequency across the 10 imputation replicates and correlate to the TRS

max_min = pd.DataFrame()
for rep in range(1,11):
    filename = 'SRTR_AA_MM_9loc_matrix_' + str(rep) + '.txt.gz'
    matrix_file = pd.read_csv(filename, sep='\t', header=0, compression='gzip')

    only_aamm = matrix_file[matrix_file.columns[(pd.Series(matrix_file.columns).str.startswith('MM'))]]

    del matrix_file  # deletes the variable from memory so only only_aamm is left since we don't need the other columns

    # Take out total MM count columns: MM_A_COUNT, etcs.
    total_MM = only_aamm[only_aamm.columns[(pd.Series(only_aamm.columns).str.endswith('COUNT'))]]
    only_aamm = only_aamm.drop(total_MM.columns, axis=1)

    mm_frequency = only_aamm.sum()  # Gives the sum of each column (each AAMM by position)

    # Sort max-min values
    max_mm = mm_frequency.sort_values(ascending=False).reset_index()

    # Rename the column headers so you can put them together into one Data Frame
    max_mm = max_mm.rename(columns={'index': 'AA_Pos_' + str(rep), 0: 'MM_Count_' + str(rep)})

    max_min = pd.concat([max_min, max_mm], axis=1)
    print(max_min.head())

# Turn into CSV for comparison
max_min.to_csv('test_aamm_frequency_imp_reps.csv', index=False, header=True)

max_min = pd.read_csv('test_aamm_frequency_imp_reps.csv', header=0)

boolean_col = max_min['AA_Pos_1'] == max_min['AA_Pos_2']
print(boolean_col.sum())
boolean_col = max_min['AA_Pos_1'] == max_min['AA_Pos_3']
print(boolean_col.sum())
boolean_col = max_min['AA_Pos_1'] == max_min['AA_Pos_4']
print(boolean_col.sum())
