#!/usr/bin/env python

import pandas as pd
from collections import defaultdict
import gzip
import matplotlib.pyplot as plt
import seaborn as sns


# Count the MMs and find the frequencies for overall MM, 1 MM count, and 2 MM count
def aamm_freq_calc(runmatch, loci_ls):
    # Dictionaries of counts and frequencies
    overallcount = defaultdict(dict)
    overallfreq = defaultdict(dict)
    onemmcount = defaultdict(dict)
    onemmfreq = defaultdict(dict)
    twommcount = defaultdict(dict)
    twommfreq = defaultdict(dict)
    for locus in loci_ls:

        aa_locus = runmatch[(runmatch['Locus'] == locus)]  # Filter to only the locus you want

        for pos in range(1, aa_locus['Position'].max() + 1):  # Look at one position at a time and count the mismatches
            aa_pos = aa_locus[aa_locus['Position'] == pos]

            mm1_count = aa_pos['1MM'].sum()
            mm2_count = aa_pos['2MM'].sum()

            if mm1_count == 0 and mm2_count == 0:
                mm_count = 0
            else:
                mm_count = (mm1_count + mm2_count) / len(aa_pos)

            if mm1_count == 0:
                mm1_freq = 0  # division by zero gives an error so have to do if statements
            else:
                mm1_freq = mm1_count / len(aa_pos)
            if mm2_count == 0:
                mm2_freq = 0
            else:
                mm2_freq = mm2_count / len(aa_pos)

            # Overall MM frequency and count
            overallcount[locus][pos] = mm1_count + mm2_count
            # print(locus + ' ' + str(pos) + ' Overall MM Count: ' + str(mm1_count + mm2_count))
            overallfreq[locus][pos] = mm_count

            # 1 MM frequency and count
            onemmcount[locus][pos] = mm1_count
            onemmfreq[locus][pos] = mm1_freq

            # 2 MM frequency and count
            twommcount[locus][pos] = mm2_count
            twommfreq[locus][pos] = mm2_freq
    return overallcount, overallfreq, onemmcount, onemmfreq, twommcount, twommfreq


# Dictionary to DataFrame function
def dict_to_df(dictionary, column_names):
    df = pd.DataFrame()
    for key, values in dictionary.items():
        for item1, item2 in values.items():
            line = pd.DataFrame({'Key': key, 'item1': item1, 'item2': item2}, index=[0])
            df = pd.concat([df, line])
    df = df.rename(columns=column_names).reset_index(drop=True)
    return df


# Turn dictionaries into DataFrames and concat the counts dataframes with the frequency dataframes
def counts_freqs(count_dict, freq_dict, which_count):
    count_df = dict_to_df(count_dict, {'Key': 'Locus', 'item1': 'Position', 'item2': which_count + '_Count'})
    freq_df = dict_to_df(freq_dict, {'Key': 'Locus', 'item1': 'Position', 'item2': which_count + '_Frequency'})

    count_freq_df = pd.merge(count_df, freq_df, how='outer',  on=['Locus', 'Position'])

    return count_freq_df


oldmatch_filename = './kidney-outcomes-sfvt/out.runmatchMC.1.txt.gz'  # Found in /lustre/project/lgragert/kidney-outcomes-sfvt
# newmatch_filename = './srtr-impute-pubsaf2306/out.runmatchMC.1.txt.gz' # '/lustre/project/lgragert/srtr-impute-pubsaf2306/out.runmatchMC.1.txt.gz'

old_runmatch = pd.read_csv(oldmatch_filename, sep='|', compression='gzip', header=None, names=['PXID', 'Locus', 'Position', '0MM', '1MM', '2MM'])
# new_runmatch = pd.read_csv(newmatch_filename, sep='|', compression='gzip', header=None, names=['PXID', 'Locus', 'Position', '0MM', '1MM', '2MM'])

# new_loci = ['A', 'C', 'B', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
old_loci = ['A', 'C', 'B', 'DRB1', 'DQB1']

# Use the def to get the MM frequency calculations for the old and new data
oldcount, oldfreq, oldcount1, oldfreq1, oldcount2, oldfreq2 = aamm_freq_calc(old_runmatch, old_loci)
del old_runmatch
# newcount, newfreq, newcount1, newfreq1, newcount2, newfreq2 = aamm_freq_calc(new_runmatch, new_loci)
# del new_runmatch  # This clears up the memory, you no longer need this DataFrame


# Turn the dictionary into a DataFrame with the frequency and counts in the same one
old_freqmm = counts_freqs(oldcount, oldfreq, 'Overall')
old_freq1mm = counts_freqs(oldcount1, oldfreq1, '1MM')
old_freq2mm = counts_freqs(oldcount2, oldfreq2, '2MM')

# new_freqmm = counts_freqs(newcount, newfreq, 'Overall')
# new_freq1mm = counts_freqs(newcount1, newfreq1, '1MM')
# new_freq2mm = counts_freqs(newcount2, newfreq2, '2MM')

# Concat all the MM counts and put into a CSV
old_freq = pd.concat([old_freqmm, old_freq1mm['1MM_Count'], old_freq1mm['1MM_Frequency'], old_freq2mm['2MM_Count'], old_freq2mm['2MM_Frequency']], axis=1)
# new_freq = pd.concat([new_freqmm, new_freq1mm['1MM_Count'], new_freq1mm['1MM_Frequency'], new_freq2mm['2MM_Count'], new_freq2mm['2MM_Frequency']], axis=1)
old_freq.to_csv('old_AAMM_frequency.csv', index=False, header=True)
# new_freq.to_csv('new_AAMM_frequency.csv', index=False, header=True)

old_freq = pd.read_csv('old_AAMM_frequency.csv', header=0)
new_freq = pd.read_csv('new_AAMM_frequency.csv', header=0)

new_freq['MM_Freq_PUBSAF2306'] = new_freq['Overall_Frequency']
old_freq['MM_Freq_PUBSAF1812'] = old_freq['Overall_Frequency']

new_freq = new_freq[['Locus', 'Position', 'MM_Freq_PUBSAF2306']]
old_freq = old_freq[['Locus', 'Position', 'MM_Freq_PUBSAF1812']]

# Filter to only old_loci for new_freq
new_old_freq = pd.merge(new_freq, old_freq, how='inner', on=['Locus', 'Position'])
new_old_freq['Abs_Diff'] = new_old_freq['MM_Freq_PUBSAF2306'] - new_old_freq['MM_Freq_PUBSAF1812']
new_old_freq.to_csv('pubsaf_absdiff_old_loci.csv', index=False, header=True)

# Keep all rows contained to show big difference in analysis
old_new_freq = pd.merge(new_freq, old_freq, how='outer', on=['Locus', 'Position'])
old_new_freq = old_new_freq.fillna(0)
old_new_freq['Abs_Diff'] = old_new_freq['MM_Freq_PUBSAF2306'] - old_new_freq['MM_Freq_PUBSAF1812']
old_new_freq.to_csv('pubsaf_absdiff_new_loci.csv', index=False, header=True)
