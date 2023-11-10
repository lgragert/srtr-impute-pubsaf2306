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

SAS files were pre-converted to tab-delimited text files by Nick Brown using R Haven.

An example R Haven SAS extract script is here `haven_sas_experiments.R` - tested on pubsaf1812

Previous version for working with pubsaf1812 extracts are located here:

https://github.com/lgragert/srtr-impute

Includes many .sas extract scripts that are no longer needed.
