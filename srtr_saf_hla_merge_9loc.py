#!/usr/bin/env python
import pandas
import numpy as np


# input file tab-delimited CSVs - converted from SRTR SAF SAS datasets using R Haven
# Using pubsaf2306

# TX_KI file contains transplant outcome info - Has the PX_ID, DONOR_ID, TX_ID, and TRR_ID
# DONOR_DECEASED file contains other locus typing - link via DON_ID
# DONOR_LIVE file contains other locus typing for living donors - link via DON_ID
# REC_HISTO file contains other locus typing - link via REC_HISTO_TX_ID

# Identifiers

# PX_ID - "Patient Identifier” - Links donor and recipient combo
# PERS_ID - "Unique Person ID for patient. Based on matches in similarity of SSN, DOB, Names and Nicknames, Gender, etc.” - Why is this different than PX_ID??? - Used as recipient identifier in ID_MAP file
# TX_ID - "Unique identifier for TX - set same for multi_org TXs"
# DON_ID - "Unique Donor ID (all donors)” - Called “DONOR_ID” in TX_KI file
# TRR_ID - "Transplant Recipient Registration (TRR) form" - "Unique identifier for TRR - unique key"
# REC_HISTO_TX_ID - "Unique identifier for Transplant - foreign key"

# load tab-delimited CSVs to dataframe 
# Engine = python prevents unicode decode errors - https://stackoverflow.com/questions/12468179/unicodedecodeerror-utf8-codec-cant-decode-byte-0x9c
tx_ki_filename = "TX_KI_decoded.txt"
# tx_ki = pandas.read_csv(tx_ki_filename, index_col="PX_ID", engine='python')
# tx_ki = pandas.read_csv(tx_ki_filename, index_col=None, engine='python')
tx_ki = pandas.read_csv(tx_ki_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t')

print ("TX_KI Loaded: " + str(len(tx_ki)))
print(list(tx_ki.columns))

donor_deceased_filename = "DONOR_DECEASED_decoded.txt"
# donor_deceased = pandas.read_csv(donor_deceased_filename, index_col="PERS_ID", engine='python')
# donor_deceased = pandas.read_csv(donor_deceased_filename, index_col=None, engine='python')
donor_deceased = pandas.read_csv(donor_deceased_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t')

print ("DONOR_DECEASED Loaded: " + str(len(donor_deceased)))
print(list(donor_deceased.columns))

rec_histo_filename = "rec_histo_decoded.txt"
# rec_histo = pandas.read_csv(rec_histo_filename, index_col="REC_HISTO_TX_ID", engine='python')
# rec_histo = pandas.read_csv(rec_histo_filename, index_col=None, engine='python')
rec_histo = pandas.read_csv(rec_histo_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t')

print ("REC_HISTO Loaded: " + str(len(rec_histo)))
print(list(rec_histo.columns))

DQA1_DPA1_HLA_filename = "SAF_DPA_DQA_decoded.txt"
DQA1_DPA1_HLA = pandas.read_csv(DQA1_DPA1_HLA_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t')

print ("DQA1 DPA1 Loaded: " + str(len(DQA1_DPA1_HLA)))
print(list(DQA1_DPA1_HLA.columns))

# subset TX_KI columns for merge
print ("Subset REC_HISTO columns: ")
tx_ki = tx_ki[['ORG_TY','PERS_ID','PX_ID','REC_TX_DT', 'REC_HISTO_TX_ID', 'REC_TX_TY',\
		'DON_TY','DON_RACE','DON_RACE_SRTR','DON_ETHNICITY_SRTR','DON_A1','DON_A2','DON_B1','DON_B2','DON_DR1','DON_DR2', \
		'REC_AGE_IN_MONTHS_AT_TX','CAN_RACE','CAN_RACE_SRTR','CAN_ETHNICITY_SRTR','REC_A1','REC_A2','REC_B1','REC_B2','REC_DR1','REC_DR2','DONOR_ID']]


# SAS filters on TX_KI
# where (REC_AGE_AT_TX ge 18) and ((DON)TY = 'C') and (ORG_TY = 'KI') and (REC_MULTI_ORG='') and (2000 le YER(REC_TX_DT) le 2017);
# drop rows that don't meet criteria from pullIDs.sas criteria

# Substitute for REC_MULTI_ORG - REC_TX_TY - "1: Single donor, single organ type TX"
# 93 REC_MULTI_ORG Char 15 $MULTORG. Text listing organs in a multi-organ TX  in TXF_HR
# 92 REC_KP_MULTI Num 8 KP_MULT. Flag indicating if a Flag indicating if a KP TX or a KP involved in a Multi

# handle dates

# convert to DateTime - already in ISO8601 format in pubsaf2306!
# tx_ki['REC_TX_DT'] = pandas.to_datetime(tx_ki['REC_TX_DT'],format='%m/%d/%y')
# filter based on dates - TODO - migrate this downstream after imputation
# tx_ki = tx_ki[(tx_ki['REC_TX_DT'] > '2000-01-01') & (tx_ki['REC_TX_DT'] < '2017-12-31')]


# Use donor and recipient race and ethnic info to NMDP Rollup Race codes

# RACE_SRTR and ETHNICITY_SRTR splits the concepts of race and ethnicity - PREFERRED
# DON_RACE and CAN_RACE is more similar to NMDP rollup race, but has nonstandard formatting for multiracial

# CAN_RACE categories
# Multi-Racial
# 8: White
# 2000: Hispanic/Latino
# 64: Asian

# CAN_RACE_SRTR categories
# WHITE: White
# MULTI: Multiracial

# CAN_RACE_ETHNICITY categories
# LATINO: Latino
# NLATIN: Non-Latino or unknown

# Using MLT category for imputation using NEMO 9-locus pipeline - HapLogic uses CAU - Consider using US

tx_ki.DON_RACE_SRTR = tx_ki.DON_RACE_SRTR.str.split(': ',expand=True)[1]
tx_ki.DON_ETHNICITY_SRTR = tx_ki.DON_ETHNICITY_SRTR.str.split(': ',expand=True)[0]
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'White') & (tx_ki.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "CAU"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Asian') & (tx_ki.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "API"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Black') & (tx_ki.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "AFA"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Native American') & (tx_ki.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "NAM"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Pacific Islander') & (tx_ki.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "API"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'White') & (tx_ki.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "HIS"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Asian') & (tx_ki.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "MLT" # MLT
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Black') & (tx_ki.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "MLT" # MLT
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Native American') & (tx_ki.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "HIS"
tx_ki.loc[((tx_ki.DON_RACE_SRTR == 'Pacific Islander') & (tx_ki.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "MLT" # MLT
tx_ki.loc[tx_ki.DON_RACE_SRTR == 'Multiracial', 'DON_RACE'] = "MLT" # MLT

tx_ki.CAN_RACE_SRTR = tx_ki.CAN_RACE_SRTR.str.split(': ',expand=True)[1]
tx_ki.CAN_ETHNICITY_SRTR = tx_ki.CAN_ETHNICITY_SRTR.str.split(': ',expand=True)[0]
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'White') & (tx_ki.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "CAU"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Asian') & (tx_ki.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "API"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Black') & (tx_ki.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "AFA"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Native American') & (tx_ki.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "NAM"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Pacific Islander') & (tx_ki.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "API"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'White') & (tx_ki.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "HIS"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Asian') & (tx_ki.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "MLT" # MLT
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Black') & (tx_ki.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "MLT" # MLT
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Native American') & (tx_ki.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "HIS"
tx_ki.loc[((tx_ki.CAN_RACE_SRTR == 'Pacific Islander') & (tx_ki.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "MLT" # MLT
tx_ki.loc[tx_ki.CAN_RACE_SRTR == 'Multiracial', 'CAN_RACE'] = "MLT" # MLT

# subset DONOR_DECEASED columns
# PERS_ID is not the same between DONOR_DECEASED and TX_KI!!
# PERS_ID range in DONOR_DECEASED is 8000000 to 8271341
# PERS_ID range in TX_KI is 2500000 to 5000000
donor_deceased = donor_deceased[['DONOR_ID','DON_C1','DON_C2','DON_DQ1','DON_DQ2','DON_DP1','DON_DP2','DON_DR51','DON_DR52','DON_DR53']]

print ("Subset DONOR_DECEASED columns: ")
print(list(donor_deceased.columns))

# subset REC_HISTO columns 
# donor typing is re-typed HLA
# rec_histo = rec_histo[['DON_A1','DON_A2','DON_B1','DON_B2','DON_CW1','DON_CW2','DON_DPW1','DON_DPW2','DON_DQW1','DON_DQW2','DON_DR1','DON_DR2','DON_DRW51','DON_DRW52','DON_DRW53', \
# 					'REC_A1','REC_A2','REC_B1','REC_B2','REC_CW1','REC_CW2','REC_DPW1','REC_DPW2','REC_DQW1','REC_DQW2','REC_DR1','REC_DR2','REC_DRW51','REC_DRW52','REC_DRW53']]

rec_histo = rec_histo[['REC_HISTO_TX_ID','REC_CW1','REC_CW2','REC_DQW1','REC_DQW2','REC_DPW1','REC_DPW2','REC_DRW51','REC_DRW52','REC_DRW53']]

print ("Subset REC_HISTO columns: ")
print(list(rec_histo.columns))

# DQA1_DPA1 rename and subset columns

# DON_ID is UNOS_ID
# DQA1_DPA1_HLA = DQA1_DPA1_HLA.rename(columns={"DON_ID": "DONOR_ID"}) 


DQA1_DPA1_HLA = DQA1_DPA1_HLA.rename(columns={"don_dqa1": "DON_DQA1", "don_dqa2": "DON_DQA2", "don_dpa1": "DON_DPA1", "don_dpa2": "DON_DPA2",
                    "cand_dqa1": "REC_DQA1", "cand_dqa2": "REC_DQA2", "cand_dpa1": "REC_DPA1", "cand_dpa2": "REC_DPA2"})
DQA1_DPA1_HLA = DQA1_DPA1_HLA[['PX_ID','DON_DQA1','DON_DQA2','DON_DPA1','DON_DPA2','REC_DQA1','REC_DQA2','REC_DPA1','REC_DPA2']]
print ("Subset DQA1_DPA1 columns: ")
print(list(DQA1_DPA1_HLA.columns))
# print (DQA1_DPA1_HLA.dtypes)
# DQA1_DPA1_HLA = DQA1_DPA1_HLA["DONOR_ID"].astype(np.int64)


# MERGE AND FILTER TABLES

print ("Merging HLA for All Transplants in SRTR SAF TX_KI table pubsaf2306: " + str(len(tx_ki)))


# select organ type kidney - excludes KP - didn't change count so commenting out
# tx_ki.loc[tx_ki.ORG_TY == 'KI: Kidney', 'ORG_TY'] = "KI"

# select only single donor, single organ type TX for REC_TX_TY
# 	Transplant Type, number of donors & organ types involved in TX
# 1: Single donor, single organ type TX
# 2: Single donor, multiple organ types TX
# 3: Multiple donors, single organ type TX  - can't use because merging DQA1 and DPA1 on patient ID
# 4: Multiple donors, multiple organ type TX
tx_ki.REC_TX_TY = tx_ki.REC_TX_TY.str.split(': ',expand=True)[1]
tx_ki.loc[(tx_ki.REC_TX_TY == "Single donor, single organ type TX"), 'REC_TX_TY'] = "1"
tx_ki = tx_ki[(tx_ki['REC_TX_TY'] == "1")]

print ("Subset to Single donor, single organ type TX: " + str(len(tx_ki)))

# merge DONOR_DECEASED on DONOR_ID to add C, DQ, and DP
tx_ki_donor_hla = tx_ki.merge(donor_deceased,how="inner",on="DONOR_ID")
print ("Merged with DONOR_DECEASED on DONOR_ID: " + str(len(tx_ki_donor_hla)))

# merge REC_HISTO on REC_HISTO_TX_ID to add C, DQ, and DP
tx_ki_donor_rec_hla = tx_ki_donor_hla.merge(rec_histo,how="inner",on="REC_HISTO_TX_ID")
print ("Merged with REC_HISTO on REC_HISTO_TX_ID: " + str(len(tx_ki_donor_rec_hla)))

# find row dropped in merge with REC_HISTO
# print (tx_ki_donor_hla[~tx_ki_donor_hla['REC_HISTO_TX_ID'].isin(rec_histo['REC_HISTO_TX_ID'])])
# print ("Missing row from REC_HISTO merge with TX_KI: ")
# print(tx_ki.loc[tx_ki['PERS_ID'] == 2679941])
# tx_ki_missing_from_rec_histo = tx_ki.loc[tx_ki['PERS_ID'] == 2679941]
# print ("Missing REC_HISTO_TX_ID: ")
# print (tx_ki_missing_from_rec_histo['REC_HISTO_TX_ID'])
# PX_ID -1656071
# PERS_ID 2679941
# REC_HISTO_TX_ID 1776506
# print (rec_histo[~rec_histo['REC_HISTO_TX_ID'].isin(rec_histo['REC_HISTO_TX_ID'])])

# merge DQA1_DPA1_HLA on PX_ID to add C, DQ, and DP
# can't merge on DONOR_ID - uses actual UNOS donor IDs which aren't in SAF
# can't use 
tx_ki_all_hla = tx_ki_donor_rec_hla.merge(DQA1_DPA1_HLA,how="left",on="PX_ID")
print ("Merged with DQA1_DPA1_HLA on PX_ID: " + str(len(tx_ki_all_hla)))


# list columns in final table
# print(list(tx_ki_all_hla.columns))

# print (tx_ki_all_hla['DON_A1'])

hla_filename = "tx_ki_hla_9loc.csv"
tx_ki_all_hla.to_csv(hla_filename, header=True, index=False)

# after SAS formats decoding, HLA field look like "203: 0203"
# these commands convert the data to what appears after the colon - "0203"
# not needed for DQA1 and DPA1 which were extracted outside of SAS
# DQA1 and DPA1 already cleanly decoded
tx_ki_all_hla.DON_A1 = tx_ki_all_hla.DON_A1.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_A2 = tx_ki_all_hla.DON_A2.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_B1 = tx_ki_all_hla.DON_B1.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_B2 = tx_ki_all_hla.DON_B2.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_C1 = tx_ki_all_hla.DON_C1.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_C2 = tx_ki_all_hla.DON_C2.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DR1 = tx_ki_all_hla.DON_DR1.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DR2 = tx_ki_all_hla.DON_DR2.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DQ1 = tx_ki_all_hla.DON_DQ1.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DQ2 = tx_ki_all_hla.DON_DQ2.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DP1 = tx_ki_all_hla.DON_DP1.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DP2 = tx_ki_all_hla.DON_DP2.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DR51 = tx_ki_all_hla.DON_DR51.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DR52 = tx_ki_all_hla.DON_DR52.str.split(': ',expand=True)[1]
tx_ki_all_hla.DON_DR53 = tx_ki_all_hla.DON_DR53.str.split(': ',expand=True)[1]

tx_ki_all_hla.REC_A1 = tx_ki_all_hla.REC_A1.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_A2 = tx_ki_all_hla.REC_A2.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_B1 = tx_ki_all_hla.REC_B1.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_B2 = tx_ki_all_hla.REC_B2.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_CW1 = tx_ki_all_hla.REC_CW1.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_CW2 = tx_ki_all_hla.REC_CW2.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DR1 = tx_ki_all_hla.REC_DR1.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DR2 = tx_ki_all_hla.REC_DR2.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DQW1 = tx_ki_all_hla.REC_DQW1.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DQW2 = tx_ki_all_hla.REC_DQW2.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DPW1 = tx_ki_all_hla.REC_DPW1.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DPW2 = tx_ki_all_hla.REC_DPW2.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DRW51 = tx_ki_all_hla.REC_DRW51.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DRW52 = tx_ki_all_hla.REC_DRW52.str.split(': ',expand=True)[1]
tx_ki_all_hla.REC_DRW53 = tx_ki_all_hla.REC_DRW53.str.split(': ',expand=True)[1]


# skip cases where no typing at A, B, or DR
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_A1'] != "Unknown"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_B1'] != "Unknown"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_DR1'] != "Unknown"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_A1'] != "Unknown"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_B1'] != "Unknown"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_DR1'] != "Unknown"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_A1'] != "Not Tested"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_B1'] != "Not Tested"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_DR1'] != "Not Tested"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_A1'] != "Not Tested"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_B1'] != "Not Tested"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_DR1'] != "Not Tested"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_A1'] != "0"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_B1'] != "0"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['DON_DR1'] != "0"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_A1'] != "0"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_B1'] != "0"]
tx_ki_all_hla = tx_ki_all_hla[tx_ki_all_hla['REC_DR1'] != "0"]

print ("Minimum Typing of A, B, DR antigens: " + str(len(tx_ki_all_hla)))

# missing DQ and DP - just skip the locus
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1 == 'Unknown', 'DON_DQ1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == 'Unknown', 'DON_DP1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1 == '0', 'DON_DQ1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == '0', 'DON_DP1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA1 == 'Unknown', 'DON_DQA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA1 == 'Unknown', 'DON_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA1 == '0', 'DON_DQA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA1 == '0', 'DON_DPA1'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW1 == 'Unknown', 'REC_DQW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == 'Unknown', 'REC_DPW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW1 == '0', 'REC_DQW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == '0', 'REC_DPW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA1 == 'Unknown', 'REC_DQA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA1 == 'Unknown', 'REC_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA1 == '0', 'REC_DQA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA1 == '0', 'REC_DPA1'] = ""

# manage homozygosity

tx_ki_all_hla.loc[tx_ki_all_hla.DON_A2 == 'No second antigen detected', 'DON_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_B2 == 'No second antigen detected', 'DON_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DR2 == 'No second antigen detected', 'DON_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA2 == 'No second antigen detected', 'DON_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ2 == 'No second antigen detected', 'DON_DQ2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == 'No second antigen detected', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == 'No second antigen detected', 'DON_DP2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_A2 == 'Not Tested', 'DON_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_B2 == 'Not Tested', 'DON_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DR2 == 'Not Tested', 'DON_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA2 == 'Not Tested', 'DON_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ2 == 'Not Tested', 'DON_DQ2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == 'Not Tested', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == 'Not Tested', 'DON_DP2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_A2 == 'Unknown', 'DON_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_B2 == 'Unknown', 'DON_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DR2 == 'Unknown', 'DON_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA2 == 'Unknown', 'DON_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ2 == 'Unknown', 'DON_DQ2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == 'Unknown', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == 'Unknown', 'DON_DP2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_A2 == '0', 'DON_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_B2 == '0', 'DON_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DR2 == '0', 'DON_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA2 == '0', 'DON_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ2 == '0', 'DON_DQ2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == '0', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == '0', 'DON_DP2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.REC_A2 == 'No second antigen detected', 'REC_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_B2 == 'No second antigen detected', 'REC_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DR2 == 'No second antigen detected', 'REC_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA2 == 'No second antigen detected', 'REC_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW2 == 'No second antigen detected', 'REC_DQW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == 'No second antigen detected', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == 'No second antigen detected', 'REC_DPW2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.REC_A2 == 'Not Tested', 'REC_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_B2 == 'Not Tested', 'REC_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DR2 == 'Not Tested', 'REC_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA2 == 'Not Tested', 'REC_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW2 == 'Not Tested', 'REC_DQW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == 'Not Tested', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == 'Not Tested', 'REC_DPW2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.REC_A2 == 'Unknown', 'REC_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_B2 == 'Unknown', 'REC_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DR2 == 'Unknown', 'REC_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA2 == 'Unknown', 'REC_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW2 == 'Unknown', 'REC_DQW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == 'Unknown', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == 'Unknown', 'REC_DPW2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.REC_A2 == '0', 'REC_A2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_B2 == '0', 'REC_B2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DR2 == '0', 'REC_DR2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA2 == '0', 'REC_DQA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW2 == '0', 'REC_DQW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == '0', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == '0', 'REC_DPW2'] = ""

# manage No 2nd antigen in first column - Some will result in two blanks but this is handled later

tx_ki_all_hla.loc[tx_ki_all_hla.DON_A1 == 'No second antigen detected', 'DON_A1'] = tx_ki_all_hla.DON_A2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_B1 == 'No second antigen detected', 'DON_B1'] = tx_ki_all_hla.DON_B2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DR1 == 'No second antigen detected', 'DON_DR1'] = tx_ki_all_hla.DON_DR2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQA1 == 'No second antigen detected', 'DON_DQA1'] = tx_ki_all_hla.DON_DQA2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1 == 'No second antigen detected', 'DON_DQ1'] = tx_ki_all_hla.DON_DQ2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA1 == 'No second antigen detected', 'DON_DPA1'] = tx_ki_all_hla.DON_DPA2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == 'No second antigen detected', 'DON_DP1'] = tx_ki_all_hla.DON_DP2

tx_ki_all_hla.loc[tx_ki_all_hla.REC_A1 == 'No second antigen detected', 'REC_A1'] = tx_ki_all_hla.REC_A2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_B1 == 'No second antigen detected', 'REC_B1'] = tx_ki_all_hla.REC_B2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DR1 == 'No second antigen detected', 'REC_DR1'] = tx_ki_all_hla.REC_DR2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQA1 == 'No second antigen detected', 'REC_DQA1'] = tx_ki_all_hla.REC_DQA2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW1 == 'No second antigen detected', 'REC_DQW1'] = tx_ki_all_hla.REC_DQW2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA1 == 'No second antigen detected', 'REC_DPA1'] = tx_ki_all_hla.REC_DPA2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == 'No second antigen detected', 'REC_DPW1'] = tx_ki_all_hla.REC_DPW2


# manage 4-digit typing without colons (DP has colons for 4-digit alelels)
 
# # https://stackoverflow.com/questions/47293729/replacing-part-of-a-string-in-pandas-column-with-and-condition
# regex=true for older Pandas version - 
tx_ki_all_hla['DON_A1'] = tx_ki_all_hla['DON_A1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_A2'] = tx_ki_all_hla['DON_A2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_B1'] = tx_ki_all_hla['DON_B1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_B2'] = tx_ki_all_hla['DON_B2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_C1'] = tx_ki_all_hla['DON_C1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_C2'] = tx_ki_all_hla['DON_C2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_DR1'] = tx_ki_all_hla['DON_DR1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_DR2'] = tx_ki_all_hla['DON_DR2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_DQ1'] = tx_ki_all_hla['DON_DQ1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['DON_DQ2'] = tx_ki_all_hla['DON_DQ2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)


tx_ki_all_hla['REC_A1'] = tx_ki_all_hla['REC_A1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_A2'] = tx_ki_all_hla['REC_A2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_B1'] = tx_ki_all_hla['REC_B1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_B2'] = tx_ki_all_hla['REC_B2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_CW1'] = tx_ki_all_hla['REC_CW1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_CW2'] = tx_ki_all_hla['REC_CW2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_DR1'] = tx_ki_all_hla['REC_DR1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_DR2'] = tx_ki_all_hla['REC_DR2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_DQW1'] = tx_ki_all_hla['REC_DQW1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_ki_all_hla['REC_DQW2'] = tx_ki_all_hla['REC_DQW2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)


# manage C serology - No second antigen can mean low expression alleles - Imputation algorithm should do this!


# manage not tested and no second antigen for C and DQ

tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == 'Not Tested', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == 'Not Tested', 'REC_CW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == 'Not Tested', 'DON_C2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == 'Not Tested', 'REC_CW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1 == 'Not Tested', 'DON_DQ1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW1 == 'Not Tested', 'REC_DQW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ2 == 'Not Tested', 'DON_DQ2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW2 == 'Not Tested', 'REC_DQW2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == 'No antigen detected', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == 'No antigen detected', 'REC_CW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1 == 'No antigen detected', 'DON_DQ1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW1 == 'No antigen detected', 'REC_DQW1'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == 'Unknown', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == 'No antigen detected', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == '0', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == 'Unknown', 'REC_CW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == 'No antigen detected', 'REC_CW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == '0', 'REC_CW1'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == 'Unknown', 'DON_C2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == 'No second antigen detected', 'DON_C2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == '0', 'DON_C2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == 'Unknown', 'REC_CW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == 'No second antigen detected', 'REC_CW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == '0', 'REC_CW2'] = ""


tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == 'No second antigen detected', 'DON_C1'] = tx_ki_all_hla.DON_C2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == 'No second antigen detected', 'REC_CW1'] = tx_ki_all_hla.REC_CW2
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1 == 'No second antigen detected', 'DON_DQ1'] = tx_ki_all_hla.DON_DQ2
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DQW1 == 'No second antigen detected', 'REC_DQW1'] = tx_ki_all_hla.REC_DQW2

# convert splits C9 and C10 to broad C3 - C9 and C10 aren't first-field DNA categories
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == '09', 'DON_C1'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == '09', 'DON_C2'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == '09', 'REC_CW1'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == '09', 'REC_CW2'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == '10', 'DON_C1'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == '10', 'DON_C2'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == '10', 'REC_CW1'] = "03"
# tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == '10', 'REC_CW2'] = "03"

# C11 - Invalid category - D20988, R4876

tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == '11', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == '11', 'DON_C2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == '11', 'REC_CW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == '11', 'REC_CW2'] = ""

# C13 - Invalid category - deleted allele

tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1 == '13', 'DON_C1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_C2 == '13', 'DON_C2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW1 == '13', 'REC_CW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_CW2 == '13', 'REC_CW2'] = ""

# convert all C typings to DNA - should already have leading zero
# tx_ki_all_hla.DON_C1 = tx_ki_all_hla['DON_C1'].str.replace(r'^(\d\d)$',r'\1:XX')
# tx_ki_all_hla.DON_C2 = tx_ki_all_hla['DON_C2'].str.replace(r'^(\d\d)$',r'\1:XX')
# tx_ki_all_hla.REC_CW1 = tx_ki_all_hla['REC_CW1'].str.replace(r'^(\d\d)$',r'\1:XX')
# tx_ki_all_hla.REC_CW2 = tx_ki_all_hla['REC_CW2'].str.replace(r'^(\d\d)$',r'\1:XX')

# this works - replaces all strings of length 4 with "FOURDIGIT"
# First part is correct - getting the substrings is the challenge
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ1.str.len() == 4, 'DON_DQ1'] = "FOURDIGIT"
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_DQ2.str.len() == 4, 'DON_DQ2'] = "FOURDIGIT"
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_A1.str.len() == 4, 'DON_A1'] = tx_ki_all_hla.loc['DON_A1'][0:1] + ":" + tx_ki_all_hla.loc['DON_A1'][2:3]
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_DR2.str.len() == 4, 'DON_DR2'] = tx_ki_all_hla['DON_DR2'][0:1] + ":" + tx_ki_all_hla['DON_DR2'][2:3]
# tx_ki_all_hla.loc[tx_ki_all_hla.DON_C1.str.len() == 4, 'DON_C1'] = "V3"


# manage DP typing - typing with "inactive" DP antigen codes (1-6) are removed
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == 'Not Tested', 'DON_DP1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == 'Not Tested', 'DON_DP2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA1 == 'Not Tested', 'DON_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == 'Not Tested', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA1 == 'Not Tested', 'REC_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == 'Not Tested', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == 'Not Tested', 'REC_DPW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == 'Not Tested', 'REC_DPW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA1 == 'No second antigen detected', 'DON_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == 'No second antigen detected', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == 'No second antigen detected', 'DON_DP1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == 'No second antigen detected', 'DON_DP2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA1 == 'No second antigen detected', 'REC_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == 'No second antigen detected', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == 'No second antigen detected', 'REC_DPW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == 'No second antigen detected', 'REC_DPW2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA1 == 'No antigen detected', 'DON_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DPA2 == 'No antigen detected', 'DON_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == 'No antigen detected', 'DON_DP1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == 'No antigen detected', 'DON_DP2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA1 == 'No antigen detected', 'REC_DPA1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPA2 == 'No antigen detected', 'REC_DPA2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == 'No antigen detected', 'REC_DPW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == 'No antigen detected', 'REC_DPW2'] = ""

tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP1 == '1 (Inactive)', 'DON_DP1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.DON_DP2 == '1 (Inactive)', 'DON_DP2'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW1 == '1 (Inactive)', 'REC_DPW1'] = ""
tx_ki_all_hla.loc[tx_ki_all_hla.REC_DPW2 == '1 (Inactive)', 'REC_DPW2'] = ""


# DRB3/4/5 typing is mostly Negative/Positive/Not Tested
# will be handled downstream in conversion to GL strings

# view Pandas data types
# print (tx_ki_all_hla.info())


# output merged HLA dataset

hla_filename = "tx_ki_hla_9loc.csv"
tx_ki_all_hla.to_csv(hla_filename, header=True, index=False)

