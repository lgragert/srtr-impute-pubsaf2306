# srtr-impute-pubsaf2306
Extract HLA typing from SRTR and format data for 9-locus high resolution HLA imputation


SRTR SAF Data Dictionary:

https://www.srtr.org/requesting-srtr-data/saf-data-dictionary/


Key data tables for HLA and covariates:

```
TX_KI - kidney transplant data
DONOR_DECEASED - deceased donor data (including HLA)
DONOR_LIVE - living donor data (including HLA)
REC_HISTO - recipient HLA data (including more HLA loci than in TX_KI)
SAF_DPA_DQA - DQA1 and DPA1 HLA typing were provided separately (donor and recip)
```

Python script for merging HLA data across SAF database files:

```
python3 srtr_saf_hla_merge_9loc.py
```


Output file with merged HLA data:

Donor/recip pairs without A or B or DR are excluded.

```
tx_ki_hla_9loc.csv

Columns:
['ORG_TY', 'PERS_ID', 'PX_ID', 'REC_TX_DT', 'REC_HISTO_TX_ID', 'DON_TY', 'DON_RACE', 'DON_RACE_SRTR', 'DON_ETHNICITY_SRTR', 'DON_A1', 'DON_A2', 'DON_B1', 'DON_B2', 'DON_DR1', 'DON_DR2', 'REC_AGE_IN_MONTHS_AT_TX', 'CAN_RACE', 'CAN_RACE_SRTR', 'CAN_ETHNICITY_SRTR', 'REC_TX_TY', 'REC_A1', 'REC_A2', 'REC_B1', 'REC_B2', 'REC_DR1', 'REC_DR2', 'DONOR_ID', 'DON_C1', 'DON_C2', 'DON_DQ1', 'DON_DQ2', 'DON_DP1', 'DON_DP2', 'DON_DR51', 'DON_DR52', 'DON_DR53', 'REC_CW1', 'REC_CW2', 'REC_DQW1', 'REC_DQW2', 'REC_DPW1', 'REC_DPW2', 'REC_DRW51', 'REC_DRW52', 'REC_DRW53', 'DON_DQA1', 'DON_DQA2', 'DON_DPA1', 'DON_DPA2', 'REC_DQA1', 'REC_DQA2', 'REC_DPA1', 'REC_DPA2']
```

Python script for converting SRTR HLA data to genotype list strings:

```
python3 srtr_hla_glstring.py
```

Output files in GLID/pull format for genotype list IDs per donor/recip:

```
glid.srtr.txt
pull.srtr.txt
```

Get population-specific cohort size in pull file to show clean conversion of race/ethnicity:

```
cut -d ',' -f2 pull.srtr.txt | sort | uniq -c

155728 AFA
30994 ASN
421387 CAU
100931 HIS
4929 MLT
5219 NAM
```

Imputation to two-field WHO alleles
Requires subdirectories with two-field NMDP haplotype freqs for various locus combos
Haplotype freqs will be provided on I2C2 after DUA is in place
Shares some logic with EM haplotype frequency estimation
Files are split by OPTN population for [AFA, ASN, CAU, HIS, MLT, NAM]
HPI is a small population (Hawaiian and Pacific Islander)
Rolled into ASN for imputation purposes because reference data too small
(1) Split GLID/pull data into population-specific files
(2) Allele list reduction using greedy algorithm - eliminates alleles not found in freqs
(3) Imputation using partition-ligation in two-locus blocks

```
perl make_slurm_split_srtr_loni.pl
sbatch ./slurm/run_slurm_srtr_split_loni.sh
perl make_slurm_greedy_srtr_loni.pl
bash ./run_sbatch_greedy_srtr_loni.sh
perl make_slurm_impute_srtr_loni.pl
bash ./run_sbatch_impute_srtr_loni.sh
```

Output of imputation pipelines - list of haplotype pairs and probabilities for each ID

```
impute.srtr.*.csv.gz 
```

Amino acid mismatch assignment:

```

```

Decoding for DQA1 and DPA1 HLA data formats came from using UNOS APIs - example JSON API query:

```
python3 dqa1_dpa1_code_api.py
```

Requires the following app and UNOS API key information (not checked included in code / repo):

```
tulane_app_id.key
unos_api_public.key
unos_api_secret.key
```

SAS files were pre-converted to tab-delimited text files by Nick Brown using R Haven.

An example R Haven SAS extract script is here `haven_sas_experiments.R` - tested on pubsaf1812

Previous version for working with pubsaf1812 extracts are located here:

https://github.com/lgragert/srtr-impute

Includes many .sas extract scripts that are no longer needed.
