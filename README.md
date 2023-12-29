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
```

Data request for DRB3/4/5 with DQA1 and DPA1

```
SAF_DPA_DQA folder - DQA1 and DPA1 HLA typing (donor and recip)
DR51-52-53 folder - DRB3/4/5 with DQA1 and DPA1 (donor and recip)
upenn_dqadpadr5153_29Nov2023.csv - File with encoded typing data 
```

Antigen code decoding files (extracted by API):

```
DR51_Lookup.csv
DR52_Lookup.csv
DR53_Lookup.csv
UNOS_typing_codes_DQA.txt
UNOS_typing_codes_DPA.txt
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

Antigen mismatch assignment:

```
python3 srtr_hla_antigen_mm.py

srtr_antigen_mm_*.csv

REC_C_MM_EQUIV_CUR - Cw antigen MM
REC_DQ_MM_EQUIV_CUR - DQ antigen MM
REC_DQA1_MM_EQUIV_CUR - DQA1 First-field MM
REC_DPA1_MM_EQUIV_CUR - DPA1 First-field MM
REC_DPB1_MM_EQUIV_CUR - DPB1 allele mismatch
```

Amino acid coordinate info by locus - starts at 1:

```
aa_matching_msf_genie.py

"A" : 341,
"B" : 338,
"C" : 342,
"DRB1" : 237,
"DRB3" : 237,
"DRB4" : 237,
"DRB5" : 237,
"DRB345" : 237,			
"DQA1" : 232,
"DQB1" : 229,
"DPA1" : 229,
"DPB1" : 229
```

Amino acid mismatch assignment:
- Full protein coordinates
- 9 loci
- DRB3/4/5 rules are complex
    - DRB3/4/5 genes are treated as if they were the same gene but different alleles. Each gene has the same number of 
    - Genotypes with one missing DRB3/4/5 gene are treated as homozygotes.
    - If recipient has no DRB3/4/5 genes, all DRB3/4/5 positions are assigned as mismatches when a DRB3/4/5 gene is present in the donor
    - An alternative approach might count mismatches among up to 4 gene copies of DRB1/3/4/5 and have up to 4 mismatches per DRB gene position.
    - Another alternative approach might consider DRB3/4/5 as separate genes and any DRB3/4/5 gene present in the donor that isn't in the recipient would cause mismatches at all positions.

Running AAMM script:
- Arguments are replicate, generateRunMatchMC, generateMatrix, generateSFVT
- SFVT columns aren't included in runs

RunMatch files generated with following commands:

```
python3 aa_mm_biopython_runmatch_genie_9loc.py 1 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 2 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 3 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 4 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 5 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 6 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 7 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 8 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 9 1 0 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 10 1 0 0
```

RunMatch file format (used as input to SAS):

```
PXID|Locus|Position|0MM|1MM|2MM
886437|A|4|1|0|0
```

RunMatch file locations:

```
/project/kamoun_shared/data_shared/srtr_impute_output/
out.runmatchMC.*.txt.gz
```

Running data matrix format:

```
python3 aa_mm_biopython_runmatch_genie_9loc.py 1 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 2 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 3 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 4 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 5 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 6 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 7 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 8 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 9 0 1 0
python3 aa_mm_biopython_runmatch_genie_9loc.py 10 0 1 0
```

Output AAMM matrix files (w/o covariates):

```
SRTR_AA_MM_9loc_matrix_*.txt
```

Slurm script for generating all replicates in parallel on Tulane Cypress supercomputer:
- Currently stalls because of pyARD MAC code download

```
aa_mm_biopython_runmatch_genie_9loc_RUN.sh
```

Joining with covariates and filtering dataset - SAS code should do this independently:

- Censoring date: 2022-12-31
- Transplant year starting 2005
- Donor age >=9
- Recipient age >=18
- First transplant only

Runs all 

```
gunzip SRTR_AA_MM_9loc_matrix_*.txt.gz
python3 construct_outcomes_vars_9loc.py
```

Number of transplant pairs after filters applied:

```
Transplant Pairs with A, B, DRB1 typing for donor/recip: 359590
Recipient Age >=18: 345363
Donor Age >=9: 331745
Transplant Date starting 2005: 212575
First Transplant Only: 185693
```

Design matrix files for FIBERS input:
- Each locus+position has its own column
- 'grffail' column set to 1 if graft failed within the first year

```
SRTR_AA_MM_9loc_matrix_grffail_*.txt.gz
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
