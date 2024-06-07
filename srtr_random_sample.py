# Get 10 random ids from the impute files
import pandas as pd

pops = ['AFA', 'ASN', 'CAU', 'HIS', 'MLT', 'NAM']
srtr = pd.DataFrame()

for pop in pops:
    file_name = 'impute.srtr.' + pop + '.csv.gz'
    pop_file = pd.read_csv(file_name, sep=',', compression='gzip', header=None)

    srtr = pd.concat([srtr, pop_file])

recip_list = srtr[srtr[0].str.startswith('R')]
donor_list = srtr[srtr[0].str.startswith('D')]

nodupes_recip = recip_list[0].drop_duplicates()
nodupes_donor = donor_list[0].drop_duplicates()

nodupes_recip = nodupes_recip.sample(n=100, random_state=6042024)
split_recip = nodupes_recip.str.split('R', expand=True)
split_recip = split_recip[1].str.split('-', expand=True)

only_recips = pd.concat([split_recip[0], split_recip[1]])
only_recips = only_recips[only_recips.str.strip().astype(bool)].reset_index(drop=True)

recips_10 = recip_list[recip_list[0].isin(nodupes_recip)]

donors_10 = pd.DataFrame()
for id in range(0,100):
    one_donor = donor_list[donor_list[0].str.endswith(only_recips.iloc[id])]
    donors_10 = pd.concat([donors_10, one_donor])

small_samp = pd.concat([recips_10, donors_10])

small_samp.to_csv('impute.srtr.AFA.100.csv.gz', index=False, header=False)
