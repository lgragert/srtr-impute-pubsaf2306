
import pandas as pd
from collections import defaultdict

# See the difference across replicates by taking the average and standard deviation of the AAMM frequency


# runmatch files only frequency calculation
def aamm_freq_runmatch(runmatch_df, loci_ls):
    # Dictionaries of counts and frequencies
    overallcount = defaultdict(dict)
    overallfreq = defaultdict(dict)
    for locus in loci_ls:

        aa_locus = runmatch_df[(runmatch_df['Locus'] == locus)]  # Filter to only the locus you want

        for pos in range(1, aa_locus['Position'].max() + 1):  # Look at one position at a time and count the mismatches
            aa_pos = aa_locus[aa_locus['Position'] == pos]

            mm1_count = aa_pos['1MM'].sum()
            mm2_count = aa_pos['2MM'].sum()

            if mm1_count == 0 and mm2_count == 0:
                mm_count = 0
            else:
                mm_count = (mm1_count + mm2_count) / len(aa_pos)

            # Overall MM frequency and count
            overallcount[locus][pos] = mm1_count + mm2_count
            # print(locus + ' ' + str(pos) + ' Overall MM Count: ' + str(mm1_count + mm2_count))
            overallfreq[locus][pos] = mm_count

    return overallcount, overallfreq


# Dictionary to DataFrame function, dictionary = {Locus: Position: MM_count}
def dict_to_df(dictionary, column_names):
    df = pd.DataFrame()
    for key, values in dictionary.items():
        for item1, item2 in values.items():
            line = pd.DataFrame({'Key': key, 'item1': item1, 'item2': item2}, index=[0])
            df = pd.concat([df, line])
    df = df.rename(columns=column_names).reset_index(drop=True)
    return df


# Turn dictionaries into DataFrames and concat the counts dataframes with the frequency dataframes
def counts_freqs(count_dict, freq_dict, which_rep):
    count_df = dict_to_df(count_dict, {'Key': 'Locus', 'item1': 'Position', 'item2': 'Overall_' + str(which_rep) + '_Count'})
    freq_df = dict_to_df(freq_dict, {'Key': 'Locus', 'item1': 'Position', 'item2': 'Overall_' + str(which_rep) + '_Frequency'})

    count_freq_df = pd.merge(count_df, freq_df, how='outer',  on=['Locus', 'Position'])

    return count_freq_df


# Calculates the average and standard deviation across replicates
def avg_stdv_calc(matrix_runmatch):
    # Separate frequency from counts to get overall and standard deviation
    overall_count = matrix_runmatch[matrix_runmatch.columns[(pd.Series(matrix_runmatch.columns).str.endswith('Count'))]]
    overall_freq = matrix_runmatch[matrix_runmatch.columns[(pd.Series(matrix_runmatch.columns).str.endswith('Frequency'))]]

    average_count = overall_count.sum(axis=1) / 10  # Since there are 10 replicates, we are averaging by that number now
    average_freq = overall_freq.sum(axis=1) / 10

    stdev_count = overall_count.std(axis=1)
    stdev_freq = overall_freq.std(axis=1)

    avg_stdv = pd.concat([matrix_runmatch['Locus'], matrix_runmatch['Position'], average_count, stdev_count,
                          average_freq, stdev_freq], axis=1)
    avg_stdv = avg_stdv.rename(columns={0: 'Average_Count', 1: 'STDV_Count', 2: "Average_Frequency", 3: 'STDV_Frequency'})

    return avg_stdv


# Matrix files
all_matrix_reps = pd.DataFrame(columns=['Locus', 'Position'])
for rep in range(1,11):
    matrix_filename = 'SRTR_AA_MM_9loc_matrix_' + str(rep) + '.txt.gz'
    matrix = pd.read_csv(matrix_filename, sep='\t', compression='gzip', header=0)
    aamm_cols = matrix[matrix.columns[(pd.Series(matrix.columns).str.startswith('MM'))]]
    del matrix
    transpose_aamm = aamm_cols.transpose()
    aamm_only = transpose_aamm[:-9]  # Last 9 rows are only *_MM_COUNTs so we do not want them
    aamm_only = aamm_only.replace(2, 1)  # so we only count all MMs as one, since we do that for runmatchMC
    aamm_only = aamm_only.sum(axis=1).reset_index()
    aamm_only[['MM', 'Locus', 'Position']] = aamm_only['index'].str.split('_', expand=True)
    # aamm_only[['MM', 'Locus_' + str(rep), 'Position_' + str(rep)]] = aamm_only['index'].str.split('_', expand=True)
    aamm_only = aamm_only[['Locus', 'Position', 0]]
    aamm_only = aamm_only.rename(columns={0: 'Overall_' + str(rep) + '_Count'})
    aamm_only['Overall_' + str(rep) + '_Frequency'] = aamm_only['Overall_' + str(rep) + '_Count'] / len(aamm_cols)

    all_matrix_reps = pd.merge(all_matrix_reps, aamm_only, how='outer', on=['Locus', 'Position'])
    print(all_matrix_reps.head())

all_matrix_reps.to_csv('matrix_aamm_frequency.csv', index=False, header=True)

# Run all previous code in Cypress, files are too large
# Get Overall Frequency and Counts by summing up all the counts and then dividing by 10 and getting the standard deviation
all_matrix_reps = pd.read_csv('matrix_aamm_frequency.csv', header=0)
average_matrix = avg_stdv_calc(all_matrix_reps)
average_matrix.to_csv('matrix_aamm_avg_stdv_reps.csv', header=True, index=False)
# Go to R to make Manhattan plots in:

# runmatchMC files
# Tricky because you have to delete the DataFrame to save some room
loci = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

all_runamtch_reps = pd.DataFrame(columns=['Locus', 'Position'])
for rep in range(1, 11):
    runmatch_filename = 'out.runmatchMC.' + str(rep) + '.txt.gz'
    runmatch = pd.read_csv(runmatch_filename, sep='|', compression='gzip', header=None)
    runmatch = runmatch.rename(columns={0: 'PXID', 1: 'Locus', 2: 'Position', 3: '0MM', 4: '1MM', 5: '2MM'})

    # Use the def to get the MM frequency calculations for runmatch
    count_runmatch, freq_runmatch = aamm_freq_runmatch(runmatch, loci)

    # Turn the dictionary into a DataFrame with the frequency and counts in the same one
    runmatch_freqmm = counts_freqs(count_runmatch, freq_runmatch, rep)

    # Merge results to get each rep for the same AA position
    all_runamtch_reps = pd.merge(all_runamtch_reps, runmatch_freqmm, how='outer', on=['Locus', 'Position'])

all_runamtch_reps.to_csv('runmatchMC_aamm_frequency.csv', index=False, header=True)


# Run all previous code in Cypress, files are too large
# Get Overall Frequency and Counts by summing up all the counts and then dividing by 10 and getting the standard deviation
all_runamtch_reps = pd.read_csv('runmatchMC_aamm_frequency.csv', header=0)
average_runmatch = avg_stdv_calc(all_runamtch_reps)
average_runmatch.to_csv('runmatch_aamm_avg_stdv_reps.csv', header=True, index=False)


