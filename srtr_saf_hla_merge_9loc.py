#!/usr/bin/env python
import pandas
import numpy as np
import os
import sys

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

sysarg_donor = sys.argv[1]
sysarg_tx_ty = sys.argv[2]
print("DON Type:", sysarg_donor)
print("TX Type:", sysarg_tx_ty)

# Load in tab-delimited TXTs to dataframes
donor_deceased_filename = "srtr2505_saf_r/donor_deceased.txt"
donor_live_filename = "srtr2505_saf_r/donor_live.txt"
rec_histo_filename = "srtr2505_saf_r/rec_histo.txt"
# Options: TX_HL, TX_HR, TX_LU, TX_KP, TX_KI, TX_PA, TX_LI, TX_IN
tx_hl_filename = "srtr2505_saf_r/tx_hl.txt"
tx_hr_filename = "srtr2505_saf_r/tx_hr.txt"
tx_lu_filename = "srtr2505_saf_r/tx_lu.txt"
tx_kp_filename = "srtr2505_saf_r/tx_kp.txt"
tx_ki_filename = "srtr2505_saf_r/tx_ki.txt"
tx_pa_filename = "srtr2505_saf_r/tx_pa.txt"
tx_li_filename = "srtr2505_saf_r/tx_li.txt"
tx_in_filename = "srtr2505_saf_r/tx_in.txt"

# Select donor type here:
donor_type = str(sysarg_donor) # Options: deceased, living
# Choose donor filepath based on type
if donor_type == "deceased":
	donor_filename = donor_deceased_filename
elif donor_type == "living":
	donor_filename = donor_live_filename
else:
	print("Invalid donor type selected. Please choose from: deceased, living")
	# Exit the script if an invalid type is selected
	exit()
print(f"{donor_filename}")

# Select TX type here: 
tx_type = str(sysarg_tx_ty) # Options: heart_lungs, heart, lungs, kidney_pancreas, kidney, pancreas, liver, intestines
# Choose TX filepath based on type
if tx_type == "heart_lungs":
	tx_filename = tx_hl_filename
elif tx_type == "heart":
	tx_filename = tx_hr_filename
elif tx_type == "lungs":
	tx_filename = tx_lu_filename
elif tx_type == "kidney_pancreas":
	tx_filename = tx_kp_filename
elif tx_type == "kidney":
	tx_filename = tx_ki_filename
elif tx_type == "pancreas":
	tx_filename = tx_pa_filename
elif tx_type == "liver":
	tx_filename = tx_li_filename
elif tx_type == "intestines":
	tx_filename = tx_in_filename
else:
	print("Invalid TX type selected. Please choose from: heart_lungs, heart, lungs, kidney_pancreas, kidney, pancreas, liver, intestines")
	# Exit the script if an invalid type is selected
	exit()
print(f"{tx_filename}")

# load tab-delimited CSVs to dataframe 
# Engine = python prevents unicode decode errors - https://stackoverflow.com/questions/12468179/unicodedecodeerror-utf8-codec-cant-decode-byte-0x9c
# tx_ki_filename = "TX_KI_decoded.txt"
# tx_ki = pandas.read_csv(tx_ki_filename, index_col="PX_ID", engine='python')
# tx_ki = pandas.read_csv(tx_ki_filename, index_col=None, engine='python')
tx_df = pandas.read_csv(tx_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t',dtype=str)

print (f"{tx_type} Loaded: " + str(len(tx_df)))
print(list(tx_df.columns))

# donor_deceased_filename = "DONOR_DECEASED_decoded.txt"
# donor_deceased = pandas.read_csv(donor_deceased_filename, index_col="PERS_ID", engine='python')
# donor_deceased = pandas.read_csv(donor_deceased_filename, index_col=None, engine='python')
donor_df = pandas.read_csv(donor_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t',dtype=str)

print (f"{donor_type} donor Loaded: " + str(len(donor_df)))
print(list(donor_df.columns))

# rec_histo_filename = "rec_histo_decoded.txt"
# rec_histo = pandas.read_csv(rec_histo_filename, index_col="REC_HISTO_TX_ID", engine='python')
# rec_histo = pandas.read_csv(rec_histo_filename, index_col=None, engine='python')
rec_histo = pandas.read_csv(rec_histo_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t',dtype=str)

print ("REC_HISTO Loaded: " + str(len(rec_histo)))
print(list(rec_histo.columns))

# Convert column headers from lowercase to uppercase
tx_df.columns = tx_df.columns.str.upper()
print (f"{tx_type} columns converted to uppercase: ")
print(list(tx_df.columns))

donor_df.columns = donor_df.columns.str.upper()
print (f"{donor_type} donor columns converted to uppercase: ")
print(list(donor_df.columns))

rec_histo.columns = rec_histo.columns.str.upper()
print ("REC_HISTO columns converted to uppercase: ")
print(list(rec_histo.columns))

DQA1_DPA1_DRB345_HLA_filename = "upenn_dqadpadr5153_29Nov2023.csv"
DQA1_DPA1_DRB345_HLA = pandas.read_csv(DQA1_DPA1_DRB345_HLA_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str)

# Name index column as "RECORD"
DQA1_DPA1_DRB345_HLA = DQA1_DPA1_DRB345_HLA.rename(columns={DQA1_DPA1_DRB345_HLA.columns[0]: 'RECORD'})

# Function to remove quotes from column names
remove_quotes = lambda x: x.replace('"', '')
# Rename columns using the remove_quotes function
DQA1_DPA1_DRB345_HLA = DQA1_DPA1_DRB345_HLA.rename(columns=remove_quotes)


print ("DQA1 DPA1 DRB345 Loaded: " + str(len(DQA1_DPA1_DRB345_HLA)))
print(list(DQA1_DPA1_DRB345_HLA.columns))

# to be commented out?
#DQA1_DPA1_DRB345_HLA_filename = "upenn_dqadpadr5153_29Nov2023.csv"
#DQA1_DPA1_DRB345_HLA = pandas.read_csv(DQA1_DPA1_DRB345_HLA_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str)


DRB3_DECODE_filename = "DR52_Lookup.csv"
DRB3_DECODE = pandas.read_csv(DRB3_DECODE_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str)
DRB4_DECODE_filename = "DR53_Lookup.csv"
DRB4_DECODE = pandas.read_csv(DRB4_DECODE_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str)
DRB5_DECODE_filename = "DR51_Lookup.csv"
DRB5_DECODE = pandas.read_csv(DRB5_DECODE_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str)

# Rename columns using the remove_quotes function
DRB3_DECODE = DRB3_DECODE.rename(columns=remove_quotes)
DRB4_DECODE = DRB4_DECODE.rename(columns=remove_quotes)
DRB5_DECODE = DRB5_DECODE.rename(columns=remove_quotes)

# Name index column as "RECORD"
DRB3_DECODE = DRB3_DECODE.rename(columns={DRB3_DECODE.columns[0]: 'RECORD'})
DRB4_DECODE = DRB4_DECODE.rename(columns={DRB4_DECODE.columns[0]: 'RECORD'})
DRB5_DECODE = DRB5_DECODE.rename(columns={DRB5_DECODE.columns[0]: 'RECORD'})

# View Decode tables
# print ("DRB3/4/5 Decode Tables:")
# print (DRB3_DECODE)
# print (DRB4_DECODE)
# print (DRB5_DECODE)

# Create a dictionary mapping values in id to descrip values
DRB3_DECODE_dict = dict(zip(DRB3_DECODE['id'], DRB3_DECODE['descrip']))
DRB4_DECODE_dict = dict(zip(DRB4_DECODE['id'], DRB4_DECODE['descrip']))
DRB5_DECODE_dict = dict(zip(DRB5_DECODE['id'], DRB5_DECODE['descrip']))

# standardize "Not Tested" - called "NT-not tested" here
DRB3_DECODE_dict["99"] = "Not Tested"
DRB4_DECODE_dict["99"] = "Not Tested"
DRB5_DECODE_dict["99"] = "Not Tested"

# print (DRB3_DECODE_dict)
# print (DRB4_DECODE_dict)
# print (DRB5_DECODE_dict)

# Confirmed DRB3/4/5 does contain allele-specific data
# print ("Unique values for DRB3/4/5")
# print(DQA1_DPA1_DRB345_HLA['DON_DRB3_1'].unique())

# DQA1 and DPA1 decoding
DQA1_DECODE_filename = "UNOS_typing_codes_DQA.txt"
DQA1_DECODE = pandas.read_csv(DQA1_DECODE_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str, delimiter="\t")
DPA1_DECODE_filename = "UNOS_typing_codes_DPA.txt"
DPA1_DECODE = pandas.read_csv(DPA1_DECODE_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str, delimiter="\t")


print (DQA1_DECODE)
print (DPA1_DECODE)


DQA1_DECODE_dict = dict(zip(DQA1_DECODE['code'], DQA1_DECODE['allele']))
DPA1_DECODE_dict = dict(zip(DPA1_DECODE['code'], DPA1_DECODE['allele']))

print (DQA1_DECODE_dict)
print (DPA1_DECODE_dict)



# DQA1_DPA1_HLA_filename = "SAF_DPA_DQA_decoded.txt"
# DQA1_DPA1_HLA = pandas.read_csv(DQA1_DPA1_HLA_filename, index_col=None, encoding='Latin-1', low_memory=False, sep='\t')

# print ("DQA1 DPA1 Loaded: " + str(len(DQA1_DPA1_HLA)))
# print(list(DQA1_DPA1_HLA.columns))

# Load in Decode tables for HLA Table Updates - 2024-08-29 (contains new info w P groups to handle 5-digit codes)
HLAUpdates_DECODE_filename = "HLA Table Updates - 2024-08-29.xlsx"

HLAUpdates_xl = pandas.ExcelFile(HLAUpdates_DECODE_filename)
print(HLAUpdates_xl.sheet_names)

# Load all sheets into a dictionary of DataFrames
HLAUpdates_xl_dict = pandas.read_excel(HLAUpdates_DECODE_filename, sheet_name=None, index_col=None, dtype=str)
# Extract specific sheets into DataFrames
HLAUpdates_A_DECODE = HLAUpdates_xl_dict['A - lkup_alocus']
HLAUpdates_B_DECODE = HLAUpdates_xl_dict['B - lkup_blocus']
HLAUpdates_C_DECODE = HLAUpdates_xl_dict['C - lkup_cwhla']
HLAUpdates_DRB1_DECODE = HLAUpdates_xl_dict['DRB1 - lkup_drlocus']
HLAUpdates_DRB5_DECODE = HLAUpdates_xl_dict['DRB5 - lkup_dr51hla']
HLAUpdates_DRB3_DECODE = HLAUpdates_xl_dict['DRB3 - lkup_dr52hla']
HLAUpdates_DRB4_DECODE = HLAUpdates_xl_dict['DRB4 - lkup_dr53hla']
HLAUpdates_DQA1_DECODE = HLAUpdates_xl_dict['DQA1 - lkup_dqahla']
HLAUpdates_DQB1_DECODE = HLAUpdates_xl_dict['DQB1 - lkup_dqhla']
HLAUpdates_DPA1_DECODE = HLAUpdates_xl_dict['DPA1 - lkup_dpahla']
HLAUpdates_DPB1_DECODE = HLAUpdates_xl_dict['DPB1 - lkup_dphla']
print ("HLAUpdates DECODE Tables Loaded: ")

# print ("HLAUpdates DECODE Table for HLA-C locus:")
# print (HLAUpdates_C_DECODE)
# print (HLAUpdates_C_DECODE.columns)

def create_hla_upd_dict(df):

	result_dict = {}
	for index, row in df.iterrows():
		id_val = str(row['id'])
		who_code = str(row['WhoCode'])
		
		# Check if id is 5-digit and starts with "22"
		if len(id_val) == 5 and id_val.startswith('22'):
			# Trim "P" from the end of WhoCode if present
			if who_code.endswith('P'):
				who_code = who_code[:-1] 
			
		if '*' in who_code:
			who_code = who_code.split('*', 1)[-1]
		
		result_dict[id_val] = who_code
	
	return result_dict

HLAUpdates_A_DECODE_dict = create_hla_upd_dict(HLAUpdates_A_DECODE)
HLAUpdates_B_DECODE_dict = create_hla_upd_dict(HLAUpdates_B_DECODE)
HLAUpdates_C_DECODE_dict = create_hla_upd_dict(HLAUpdates_C_DECODE)
HLAUpdates_DRB1_DECODE_dict = create_hla_upd_dict(HLAUpdates_DRB1_DECODE)
HLAUpdates_DRB5_DECODE_dict = create_hla_upd_dict(HLAUpdates_DRB5_DECODE)
HLAUpdates_DRB3_DECODE_dict = create_hla_upd_dict(HLAUpdates_DRB3_DECODE)
HLAUpdates_DRB4_DECODE_dict = create_hla_upd_dict(HLAUpdates_DRB4_DECODE)
HLAUpdates_DQA1_DECODE_dict = create_hla_upd_dict(HLAUpdates_DQA1_DECODE)
HLAUpdates_DQB1_DECODE_dict = create_hla_upd_dict(HLAUpdates_DQB1_DECODE)
HLAUpdates_DPA1_DECODE_dict = create_hla_upd_dict(HLAUpdates_DPA1_DECODE)
HLAUpdates_DPB1_DECODE_dict = create_hla_upd_dict(HLAUpdates_DPB1_DECODE)

# Decoding the "Null_allele" 
HLAUpdates_A_DECODE_dict["10000"] = "Null_allele" 
HLAUpdates_B_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_C_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DRB1_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DRB5_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DRB3_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DRB4_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DQA1_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DQB1_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DPA1_DECODE_dict["10000"] = "Null_allele"
HLAUpdates_DPB1_DECODE_dict["10000"] = "Null_allele"

print ("HLAUpdates DECODE Dictionaries Created: ")
print ("HLAUpdates A DECODE Dictionary: ")
print(HLAUpdates_A_DECODE_dict)

# HLAUpdates_DECODE_dicts = {
# 	'A': HLAUpdates_A_DECODE_dict,
# 	'B': HLAUpdates_B_DECODE_dict,
# 	'C': HLAUpdates_C_DECODE_dict,
# 	'DRB1': HLAUpdates_DRB1_DECODE_dict,
# 	'DRB5': HLAUpdates_DRB5_DECODE_dict,
# 	'DRB3': HLAUpdates_DRB3_DECODE_dict,
# 	'DRB4': HLAUpdates_DRB4_DECODE_dict,
# 	'DQA1': HLAUpdates_DQA1_DECODE_dict,
# 	'DQB1': HLAUpdates_DQB1_DECODE_dict,
# 	'DPA1': HLAUpdates_DPA1_DECODE_dict,
# 	'DPB1': HLAUpdates_DPB1_DECODE_dict
# }

# for locus, decode_dict in HLAUpdates_DECODE_dicts.items():
#     filtered_dict = {k: v for k, v in decode_dict.items() if (len(k) == 5 or len(k) == 3)}
#     HLAUpdates_DECODE_dicts[locus] = filtered_dict

HLAUpdates_A_DECODE_dict = {k: v for k, v in HLAUpdates_A_DECODE_dict.items() if (len(k) == 5 or len(k) == 3) and k != '210'}
HLAUpdates_B_DECODE_dict = {k: v for k, v in HLAUpdates_B_DECODE_dict.items() if (len(k) == 5 or len(k) == 3) and k != '703'}
HLAUpdates_C_DECODE_dict = {k: v for k, v in HLAUpdates_C_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DRB1_DECODE_dict = {k: v for k, v in HLAUpdates_DRB1_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DRB5_DECODE_dict = {k: v for k, v in HLAUpdates_DRB5_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DRB3_DECODE_dict = {k: v for k, v in HLAUpdates_DRB3_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DRB4_DECODE_dict = {k: v for k, v in HLAUpdates_DRB4_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DQA1_DECODE_dict = {k: v for k, v in HLAUpdates_DQA1_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DQB1_DECODE_dict = {k: v for k, v in HLAUpdates_DQB1_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DPA1_DECODE_dict = {k: v for k, v in HLAUpdates_DPA1_DECODE_dict.items() if (len(k) == 5 or len(k) == 3)}
HLAUpdates_DPB1_DECODE_dict = {k: v for k, v in HLAUpdates_DPB1_DECODE_dict.items() if (len(k) == 5 or len(k) == 3) and k not in ['106', '102']}

print ("HLAUpdates DECODE Dictionaries Filtered: ")
print ("HLAUpdates A DECODE Dictionary: ")
print(HLAUpdates_A_DECODE_dict)
# print ("HLAUpdates B DECODE Dictionary: ")
# print(HLAUpdates_B_DECODE_dict)



# subset TX_KI columns for merge
print (f"Subset {tx_type} columns: ")
tx_df = tx_df[['ORG_TY','PERS_ID','PX_ID','REC_TX_DT', 'REC_HISTO_TX_ID', 'REC_TX_TY',\
		'DON_TY','DON_RACE','DON_RACE_SRTR','DON_ETHNICITY_SRTR','DON_A1','DON_A2','DON_B1','DON_B2','DON_DR1','DON_DR2', \
		'REC_AGE_IN_MONTHS_AT_TX','CAN_RACE','CAN_RACE_SRTR','CAN_ETHNICITY_SRTR','REC_A1','REC_A2','REC_B1','REC_B2','REC_DR1','REC_DR2','DONOR_ID']]

print (f"Subset {tx_type} columns: ")
print(list(tx_df.columns))


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

# DON_RACE and CAN_RACE categories
# Multi-Racial
# 8: White
# 16: Black or African American
# 32: American Indian or Alaska Native
# 64: Asian
# 128: Native Hawaiian or Other Pacific Islander
# 1024: Unknown (for Donor Referral only)
# 2000: Hispanic/Latino

# DON_RACE_SRTR
# ASIAN	ASIAN: Asian
# BLACK	BLACK: Black
# MULTI	MULTI: Multiracial
# NATIVE	NATIVE: Native American
# PACIFIC	PACIFIC: Pacific Islander
# WHITE	WHITE: White

# ETHNICITY categories
# LATINO: Latino
# NLATIN: Non-Latino or unknown

# Using MLT category for imputation using NEMO 9-locus pipeline
# TODO - Consider using Overall US instead of MLT to impute MLT but hasn't been run yet for two-field

# PERS_ID, PX_ID, REC_HISTO_TX_ID, DONOR_ID all loaded as int64 instead of strings
# print (tx_ki.dtypes)

# print data frame for PX_IDs to track down NAs
# result = tx_ki.loc[tx_ki.PX_ID.astype('str') == "347820"]
# print (result['DON_RACE'])
# print (result['CAN_RACE'])
# print (result['DON_RACE_SRTR'])
# print (result['CAN_RACE_SRTR'])

# abbreviate organ type data to "KI"
#tx_ki.loc[tx_ki.ORG_TY == 'KI: Kidney', 'ORG_TY'] = "KI" ### Commenting out b/c pubsaf2505 already has "KI"



#tx_ki.DON_RACE_SRTR = tx_ki.DON_RACE_SRTR.str.split(': ',expand=True)[1]
#tx_ki.DON_RACE = tx_ki.DON_RACE.str.split(': ',expand=True)[1]
#tx_ki.DON_ETHNICITY_SRTR = tx_ki.DON_ETHNICITY_SRTR.str.split(': ',expand=True)[0]
tx_df.loc[((tx_df.DON_RACE_SRTR == 'WHITE') & (tx_df.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "CAU"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'ASIAN') & (tx_df.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "ASN"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'BLACK') & (tx_df.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "AFA"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'NATIVE') & (tx_df.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "NAM"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'PACIFIC') & (tx_df.DON_ETHNICITY_SRTR == 'NLATIN')), 'DON_RACE'] = "ASN"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'WHITE') & (tx_df.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "HIS"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'ASIAN') & (tx_df.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "MLT" # MLT
tx_df.loc[((tx_df.DON_RACE_SRTR == 'BLACK') & (tx_df.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "MLT" # MLT
tx_df.loc[((tx_df.DON_RACE_SRTR == 'NATIVE') & (tx_df.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "HIS"
tx_df.loc[((tx_df.DON_RACE_SRTR == 'PACIFIC') & (tx_df.DON_ETHNICITY_SRTR == 'LATINO')), 'DON_RACE'] = "MLT" # MLT
tx_df.loc[tx_df.DON_RACE_SRTR == 'MULTI', 'DON_RACE'] = "MLT" # MLT
#tx_ki.loc[tx_ki.DON_RACE == 'Unknown (for Donor Referral only)', 'DON_RACE'] = "MLT" # MLT
tx_df.loc[tx_df.DON_RACE_SRTR.isna(), 'DON_RACE'] = "MLT" # MLT

#tx_ki.CAN_RACE_SRTR = tx_ki.CAN_RACE_SRTR.str.split(': ',expand=True)[1]
#tx_ki.CAN_ETHNICITY_SRTR = tx_ki.CAN_ETHNICITY_SRTR.str.split(': ',expand=True)[0]
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'WHITE') & (tx_df.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "CAU"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'ASIAN') & (tx_df.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "ASN"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'BLACK') & (tx_df.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "AFA"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'NATIVE') & (tx_df.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "NAM"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'PACIFIC') & (tx_df.CAN_ETHNICITY_SRTR == 'NLATIN')), 'CAN_RACE'] = "ASN"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'WHITE') & (tx_df.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "HIS"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'ASIAN') & (tx_df.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "MLT" # MLT
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'BLACK') & (tx_df.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "MLT" # MLT
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'NATIVE') & (tx_df.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "HIS"
tx_df.loc[((tx_df.CAN_RACE_SRTR == 'PACIFIC') & (tx_df.CAN_ETHNICITY_SRTR == 'LATINO')), 'CAN_RACE'] = "MLT" # MLT
tx_df.loc[tx_df.CAN_RACE_SRTR == 'MULTI', 'CAN_RACE'] = "MLT" # MLT
tx_df.loc[tx_df.CAN_RACE_SRTR.isna(), 'CAN_RACE'] = "MLT" # MLT

# subset DONOR_DECEASED columns
# PERS_ID is not the same between DONOR_DECEASED and TX_KI!!
# PERS_ID range in DONOR_DECEASED is 8000000 to 8271341
# PERS_ID range in TX_KI is 2500000 to 5000000
donor_df = donor_df[['DONOR_ID','DON_C1','DON_C2','DON_DQ1','DON_DQ2','DON_DP1','DON_DP2','DON_DR51','DON_DR52','DON_DR53']]

print (f"Subset {donor_type} donor columns: ")
print(list(donor_df.columns))

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


# DQA1_DPA1_HLA = DQA1_DPA1_HLA.rename(columns={"don_dqa1": "DON_DQA1", "don_dqa2": "DON_DQA2", "don_dpa1": "DON_DPA1", "don_dpa2": "DON_DPA2",
#                     "cand_dqa1": "REC_DQA1", "cand_dqa2": "REC_DQA2", "cand_dpa1": "REC_DPA1", "cand_dpa2": "REC_DPA2"})
# DQA1_DPA1_HLA = DQA1_DPA1_HLA[['PX_ID','DON_DQA1','DON_DQA2','DON_DPA1','DON_DPA2','REC_DQA1','REC_DQA2','REC_DPA1','REC_DPA2']]
# print ("Subset DQA1_DPA1 columns: ")
# print(list(DQA1_DPA1_HLA.columns))

print (DQA1_DPA1_DRB345_HLA.dtypes)

DQA1_DPA1_DRB345_HLA = DQA1_DPA1_DRB345_HLA.rename(
    columns={"don_dqa1": "DON_DQA1", "don_dqa2": "DON_DQA2", \
			"don_dpa1": "DON_DPA1", "don_dpa2": "DON_DPA2", \
			"cand_dqa1": "REC_DQA1", "cand_dqa2": "REC_DQA2", \
			"cand_dpa1": "REC_DPA1", "cand_dpa2": "REC_DPA2", \
			"don_dr51": "DON_DRB5_1", "don_dr51_2": "DON_DRB5_2", \
			"don_dr_52": "DON_DRB3_1", "don_dr_52_2": "DON_DRB3_2", \
			"don_dr53": "DON_DRB4_1", "don_dr53_2": "DON_DRB4_2", \
			"cand_dr51": "REC_DRB5_1", "cand_dr51_2": "REC_DRB5_2", \
			"cand_dr52": "REC_DRB3_1", "cand_dr52_2": "REC_DRB3_2", \
			"cand_dr53": "REC_DRB4_1", "cand_dr53_2": "REC_DRB4_2"})


DQA1_DPA1_DRB345_HLA = \
 	DQA1_DPA1_DRB345_HLA[['PX_ID', \
     'DON_DRB3_1','DON_DRB3_2','DON_DRB4_1','DON_DRB4_2','DON_DRB5_1','DON_DRB5_2', \
     'DON_DQA1','DON_DQA2','DON_DPA1','DON_DPA2', \
 	 'REC_DRB3_1','REC_DRB3_2','REC_DRB4_1','REC_DRB4_2','REC_DRB5_1','REC_DRB5_2', \
     'REC_DQA1','REC_DQA2','REC_DPA1','REC_DPA2' \
     ]]
print ("Subset DQA1_DPA1_DRB345 columns: ")
print(list(DQA1_DPA1_DRB345_HLA.columns))



# MERGE AND FILTER TABLES

print (f"All Transplants in SRTR SAF {tx_type} table pubsaf2505: " + str(len(tx_df)))


# select organ type kidney - excludes KP - didn't change count so commenting out
# tx_ki.loc[tx_ki.ORG_TY == 'KI: Kidney', 'ORG_TY'] = "KI"

# select only single donor, single organ type TX for REC_TX_TY
# 	Transplant Type, number of donors & organ types involved in TX
# 1: Single donor, single organ type TX
# 2: Single donor, multiple organ types TX
# 3: Multiple donors, single organ type TX  - can't use because merging DQA1 and DPA1 on patient ID
# 4: Multiple donors, multiple organ type TX

print(tx_df['REC_TX_TY'].unique())

#tx_ki.REC_TX_TY = tx_ki.REC_TX_TY.str.split(': ',expand=True)[1]
#tx_ki.loc[(tx_ki.REC_TX_TY == "Single donor, single organ type TX"), 'REC_TX_TY'] = "1"
tx_ki = tx_ki[(tx_ki['REC_TX_TY'] == "1")]

print ("Subset to Single donor, single organ type TX: " + str(len(tx_df)))

# merge DONOR_DECEASED/DONOR_LIVE on DONOR_ID to add C, DQ, and DP
tx_df_donor_hla = tx_df.merge(donor_df,how="inner",on="DONOR_ID")
print (f"Inner join with {donor_type} donor on DONOR_ID: " + str(len(tx_df_donor_hla)))

# Add row to REC_HISTO with missing REC_HISTO_TX_ID
missing_rec_histo_id = {'REC_HISTO_TX_ID': "1776506"}
rec_histo = rec_histo._append(missing_rec_histo_id, ignore_index=True) # uses private method

# merge REC_HISTO on REC_HISTO_TX_ID to add C, DQ, and DP
tx_df_donor_rec_hla = tx_df_donor_hla.merge(rec_histo,how="inner",on="REC_HISTO_TX_ID")
print ("Inner join with REC_HISTO on REC_HISTO_TX_ID: " + str(len(tx_df_donor_rec_hla)))

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
tx_df_all_hla = tx_df_donor_rec_hla.merge(DQA1_DPA1_DRB345_HLA,how="left",on="PX_ID")
print ("Left join with DQA1_DPA1_DRB345_HLA on PX_ID: " + str(len(tx_df_all_hla)))


# list columns in final table
# print(list(tx_ki_all_hla.columns))

# print (tx_ki_all_hla['DON_A1'])

hla_filename = f"{donor_type}_{tx_type}_hla_9loc.csv"
tx_df_all_hla.to_csv(hla_filename, header=True, index=False)
print(f"{donor_type} {tx_type} HLA data saved to: {hla_filename}")


# DECODE DQA1, DPA1, DRB3/4/5
tx_df_all_hla['DON_DRB3_1'] = tx_df_all_hla['DON_DRB3_1'].map(DRB3_DECODE_dict)
tx_df_all_hla['DON_DRB3_2'] = tx_df_all_hla['DON_DRB3_2'].map(DRB3_DECODE_dict)
tx_df_all_hla['DON_DRB4_1'] = tx_df_all_hla['DON_DRB4_1'].map(DRB4_DECODE_dict)
tx_df_all_hla['DON_DRB4_2'] = tx_df_all_hla['DON_DRB4_2'].map(DRB4_DECODE_dict)
tx_df_all_hla['DON_DRB5_1'] = tx_df_all_hla['DON_DRB5_1'].map(DRB5_DECODE_dict)
tx_df_all_hla['DON_DRB5_2'] = tx_df_all_hla['DON_DRB5_2'].map(DRB5_DECODE_dict)

tx_df_all_hla['REC_DRB3_1'] = tx_df_all_hla['REC_DRB3_1'].map(DRB3_DECODE_dict)
tx_df_all_hla['REC_DRB3_2'] = tx_df_all_hla['REC_DRB3_2'].map(DRB3_DECODE_dict)
tx_df_all_hla['REC_DRB4_1'] = tx_df_all_hla['REC_DRB4_1'].map(DRB4_DECODE_dict)
tx_df_all_hla['REC_DRB4_2'] = tx_df_all_hla['REC_DRB4_2'].map(DRB4_DECODE_dict)
tx_df_all_hla['REC_DRB5_1'] = tx_df_all_hla['REC_DRB5_1'].map(DRB5_DECODE_dict)
tx_df_all_hla['REC_DRB5_2'] = tx_df_all_hla['REC_DRB5_2'].map(DRB5_DECODE_dict)

print ("Unique values for DRB3/4/5")
print(tx_df_all_hla['DON_DRB3_1'].unique())
print(tx_df_all_hla['DON_DRB5_1'].unique())

print(tx_df_all_hla['DON_DQA1'].unique())

tx_df_all_hla['DON_DQA1'] = tx_df_all_hla['DON_DQA1'].map(DQA1_DECODE_dict)
tx_df_all_hla['DON_DQA2'] = tx_df_all_hla['DON_DQA2'].map(DQA1_DECODE_dict)
tx_df_all_hla['DON_DPA1'] = tx_df_all_hla['DON_DPA1'].map(DPA1_DECODE_dict)
tx_df_all_hla['DON_DPA2'] = tx_df_all_hla['DON_DPA2'].map(DPA1_DECODE_dict)

tx_df_all_hla['REC_DQA1'] = tx_df_all_hla['REC_DQA1'].map(DQA1_DECODE_dict)
tx_df_all_hla['REC_DQA2'] = tx_df_all_hla['REC_DQA2'].map(DQA1_DECODE_dict)
tx_df_all_hla['REC_DPA1'] = tx_df_all_hla['REC_DPA1'].map(DPA1_DECODE_dict)
tx_df_all_hla['REC_DPA2'] = tx_df_all_hla['REC_DPA2'].map(DPA1_DECODE_dict)

print ("Unique values for DQA1 and DPA1")
print(tx_df_all_hla['DON_DQA1'].unique())
print(tx_df_all_hla['DON_DPA1'].unique())

# %% 

# Decode HLA Updates (5-digit P codes) for 
# A, B, C, DRB1, DRB5, DRB3, DRB4, DQA1, DQB1, DPA1, DPB1 

# tx_df_all_hla['DON_A1'] = tx_df_all_hla['DON_A1'].map(HLAUpdates_A_DECODE_dict)
# tx_df_all_hla['DON_A2'] = tx_df_all_hla['DON_A2'].map(HLAUpdates_A_DECODE_dict)
# tx_df_all_hla['DON_B1'] = tx_df_all_hla['DON_B1'].map(HLAUpdates_B_DECODE_dict)
# tx_df_all_hla['DON_B2'] = tx_df_all_hla['DON_B2'].map(HLAUpdates_B_DECODE_dict)
# tx_df_all_hla['DON_C1'] = tx_df_all_hla['DON_C1'].map(HLAUpdates_C_DECODE_dict)
# tx_df_all_hla['DON_C2'] = tx_df_all_hla['DON_C2'].map(HLAUpdates_C_DECODE_dict)
# tx_df_all_hla['DON_DR1'] = tx_df_all_hla['DON_DR1'].map(HLAUpdates_DRB1_DECODE_dict)
# tx_df_all_hla['DON_DR2'] = tx_df_all_hla['DON_DR2'].map(HLAUpdates_DRB1_DECODE_dict)
# tx_df_all_hla['DON_DR51'] = tx_df_all_hla['DON_DR51'].map(HLAUpdates_DRB5_DECODE_dict) # TODO - DON_DRB5_1 and DON_DRB5_2 ?
# tx_df_all_hla['DON_DR52'] = tx_df_all_hla['DON_DR52'].map(HLAUpdates_DRB3_DECODE_dict) # TODO - DON_DRB3_1 and DON_DRB3_2 ?
# tx_df_all_hla['DON_DR53'] = tx_df_all_hla['DON_DR53'].map(HLAUpdates_DRB4_DECODE_dict) # TODO - DON_DRB4_1 and DON_DRB4_2 ?
# tx_df_all_hla['DON_DQA1'] = tx_df_all_hla['DON_DQA1'].map(HLAUpdates_DQA1_DECODE_dict)
# tx_df_all_hla['DON_DQA2'] = tx_df_all_hla['DON_DQA2'].map(HLAUpdates_DQA1_DECODE_dict)
# tx_df_all_hla['DON_DQ1'] = tx_df_all_hla['DON_DQ1'].map(HLAUpdates_DQB1_DECODE_dict)
# tx_df_all_hla['DON_DQ2'] = tx_df_all_hla['DON_DQ2'].map(HLAUpdates_DQB1_DECODE_dict)
# tx_df_all_hla['DON_DPA1'] = tx_df_all_hla['DON_DPA1'].map(HLAUpdates_DPA1_DECODE_dict)
# tx_df_all_hla['DON_DPA2'] = tx_df_all_hla['DON_DPA2'].map(HLAUpdates_DPA1_DECODE_dict)
# tx_df_all_hla['DON_DP1'] = tx_df_all_hla['DON_DP1'].map(HLAUpdates_DPB1_DECODE_dict)
# tx_df_all_hla['DON_DP2'] = tx_df_all_hla['DON_DP2'].map(HLAUpdates_DPB1_DECODE_dict)

# tx_df_all_hla['REC_A1'] = tx_df_all_hla['REC_A1'].map(HLAUpdates_A_DECODE_dict)
# tx_df_all_hla['REC_A2'] = tx_df_all_hla['REC_A2'].map(HLAUpdates_A_DECODE_dict)
# tx_df_all_hla['REC_B1'] = tx_df_all_hla['REC_B1'].map(HLAUpdates_B_DECODE_dict)
# tx_df_all_hla['REC_B2'] = tx_df_all_hla['REC_B2'].map(HLAUpdates_B_DECODE_dict)
# tx_df_all_hla['REC_CW1'] = tx_df_all_hla['REC_CW1'].map(HLAUpdates_C_DECODE_dict)
# tx_df_all_hla['REC_CW2'] = tx_df_all_hla['REC_CW2'].map(HLAUpdates_C_DECODE_dict)
# tx_df_all_hla['REC_DR1'] = tx_df_all_hla['REC_DR1'].map(HLAUpdates_DRB1_DECODE_dict)
# tx_df_all_hla['REC_DR2'] = tx_df_all_hla['REC_DR2'].map(HLAUpdates_DRB1_DECODE_dict)
# tx_df_all_hla['REC_DRW51'] = tx_df_all_hla['REC_DRW51'].map(HLAUpdates_DRB5_DECODE_dict) # TODO - REC_DRB5_1 and REC_DRB5_2 ?
# tx_df_all_hla['REC_DRW52'] = tx_df_all_hla['REC_DRW52'].map(HLAUpdates_DRB3_DECODE_dict) # TODO - REC_DRB3_1 and REC_DRB3_2 ?
# tx_df_all_hla['REC_DRW53'] = tx_df_all_hla['REC_DRW53'].map(HLAUpdates_DRB4_DECODE_dict) # TODO - REC_DRB4_1 and REC_DRB4_2 ?
# tx_df_all_hla['REC_DQA1'] = tx_df_all_hla['REC_DQA1'].map(HLAUpdates_DQA1_DECODE_dict)
# tx_df_all_hla['REC_DQA2'] = tx_df_all_hla['REC_DQA2'].map(HLAUpdates_DQA1_DECODE_dict)
# tx_df_all_hla['REC_DQW1'] = tx_df_all_hla['REC_DQW1'].map(HLAUpdates_DQB1_DECODE_dict)
# tx_df_all_hla['REC_DQW2'] = tx_df_all_hla['REC_DQW2'].map(HLAUpdates_DQB1_DECODE_dict)
# tx_df_all_hla['REC_DPA1'] = tx_df_all_hla['REC_DPA1'].map(HLAUpdates_DPA1_DECODE_dict)
# tx_df_all_hla['REC_DPA2'] = tx_df_all_hla['REC_DPA2'].map(HLAUpdates_DPA1_DECODE_dict)
# tx_df_all_hla['REC_DPW1'] = tx_df_all_hla['REC_DPW1'].map(HLAUpdates_DPB1_DECODE_dict)
# tx_df_all_hla['REC_DPW2'] = tx_df_all_hla['REC_DPW2'].map(HLAUpdates_DPB1_DECODE_dict)

tx_df_all_hla['DON_A1'] = tx_df_all_hla['DON_A1'].apply(lambda x: HLAUpdates_A_DECODE_dict.get(x, x))
tx_df_all_hla['DON_A2'] = tx_df_all_hla['DON_A2'].apply(lambda x: HLAUpdates_A_DECODE_dict.get(x, x))
tx_df_all_hla['DON_B1'] = tx_df_all_hla['DON_B1'].apply(lambda x: HLAUpdates_B_DECODE_dict.get(x, x))
tx_df_all_hla['DON_B2'] = tx_df_all_hla['DON_B2'].apply(lambda x: HLAUpdates_B_DECODE_dict.get(x, x))
tx_df_all_hla['DON_C1'] = tx_df_all_hla['DON_C1'].apply(lambda x: HLAUpdates_C_DECODE_dict.get(x, x))
tx_df_all_hla['DON_C2'] = tx_df_all_hla['DON_C2'].apply(lambda x: HLAUpdates_C_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DR1'] = tx_df_all_hla['DON_DR1'].apply(lambda x: HLAUpdates_DRB1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DR2'] = tx_df_all_hla['DON_DR2'].apply(lambda x: HLAUpdates_DRB1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DR51'] = tx_df_all_hla['DON_DR51'].apply(lambda x: HLAUpdates_DRB5_DECODE_dict.get(x, x)) # TODO - DON_DRB5_1 and DON_DRB5_2 ?
tx_df_all_hla['DON_DR52'] = tx_df_all_hla['DON_DR52'].apply(lambda x: HLAUpdates_DRB3_DECODE_dict.get(x, x)) # TODO - DON_DRB3_1 and DON_DRB3_2 ?
tx_df_all_hla['DON_DR53'] = tx_df_all_hla['DON_DR53'].apply(lambda x: HLAUpdates_DRB4_DECODE_dict.get(x, x)) # TODO - DON_DRB4_1 and DON_DRB4_2 ?
tx_df_all_hla['DON_DQA1'] = tx_df_all_hla['DON_DQA1'].apply(lambda x: HLAUpdates_DQA1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DQA2'] = tx_df_all_hla['DON_DQA2'].apply(lambda x: HLAUpdates_DQA1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DQ1'] = tx_df_all_hla['DON_DQ1'].apply(lambda x: HLAUpdates_DQB1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DQ2'] = tx_df_all_hla['DON_DQ2'].apply(lambda x: HLAUpdates_DQB1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DPA1'] = tx_df_all_hla['DON_DPA1'].apply(lambda x: HLAUpdates_DPA1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DPA2'] = tx_df_all_hla['DON_DPA2'].apply(lambda x: HLAUpdates_DPA1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DP1'] = tx_df_all_hla['DON_DP1'].apply(lambda x: HLAUpdates_DPB1_DECODE_dict.get(x, x))
tx_df_all_hla['DON_DP2'] = tx_df_all_hla['DON_DP2'].apply(lambda x: HLAUpdates_DPB1_DECODE_dict.get(x, x))

tx_df_all_hla['REC_A1'] = tx_df_all_hla['REC_A1'].apply(lambda x: HLAUpdates_A_DECODE_dict.get(x, x))
tx_df_all_hla['REC_A2'] = tx_df_all_hla['REC_A2'].apply(lambda x: HLAUpdates_A_DECODE_dict.get(x, x))
tx_df_all_hla['REC_B1'] = tx_df_all_hla['REC_B1'].apply(lambda x: HLAUpdates_B_DECODE_dict.get(x, x))
tx_df_all_hla['REC_B2'] = tx_df_all_hla['REC_B2'].apply(lambda x: HLAUpdates_B_DECODE_dict.get(x, x))
tx_df_all_hla['REC_CW1'] = tx_df_all_hla['REC_CW1'].apply(lambda x: HLAUpdates_C_DECODE_dict.get(x, x))
tx_df_all_hla['REC_CW2'] = tx_df_all_hla['REC_CW2'].apply(lambda x: HLAUpdates_C_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DR1'] = tx_df_all_hla['REC_DR1'].apply(lambda x: HLAUpdates_DRB1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DR2'] = tx_df_all_hla['REC_DR2'].apply(lambda x: HLAUpdates_DRB1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DRW51'] = tx_df_all_hla['REC_DRW51'].apply(lambda x: HLAUpdates_DRB5_DECODE_dict.get(x, x)) # TODO - REC_DRB5_1 and REC_DRB5_2 ?
tx_df_all_hla['REC_DRW52'] = tx_df_all_hla['REC_DRW52'].apply(lambda x: HLAUpdates_DRB3_DECODE_dict.get(x, x)) # TODO - REC_DRB3_1 and REC_DRB3_2 ?
tx_df_all_hla['REC_DRW53'] = tx_df_all_hla['REC_DRW53'].apply(lambda x: HLAUpdates_DRB4_DECODE_dict.get(x, x)) # TODO - REC_DRB4_1 and REC_DRB4_2 ?
tx_df_all_hla['REC_DQA1'] = tx_df_all_hla['REC_DQA1'].apply(lambda x: HLAUpdates_DQA1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DQA2'] = tx_df_all_hla['REC_DQA2'].apply(lambda x: HLAUpdates_DQA1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DQW1'] = tx_df_all_hla['REC_DQW1'].apply(lambda x: HLAUpdates_DQB1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DQW2'] = tx_df_all_hla['REC_DQW2'].apply(lambda x: HLAUpdates_DQB1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DPA1'] = tx_df_all_hla['REC_DPA1'].apply(lambda x: HLAUpdates_DPA1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DPA2'] = tx_df_all_hla['REC_DPA2'].apply(lambda x: HLAUpdates_DPA1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DPW1'] = tx_df_all_hla['REC_DPW1'].apply(lambda x: HLAUpdates_DPB1_DECODE_dict.get(x, x))
tx_df_all_hla['REC_DPW2'] = tx_df_all_hla['REC_DPW2'].apply(lambda x: HLAUpdates_DPB1_DECODE_dict.get(x, x))

# after SAS formats decoding, HLA field look like "203: 0203"
# these commands convert the data to what appears after the colon - "0203"
# not needed for DQA1 and DPA1 which were extracted outside of SAS
# DQA1 and DPA1 already cleanly decoded
# tx_ki_all_hla.DON_A1 = tx_ki_all_hla.DON_A1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_A2 = tx_ki_all_hla.DON_A2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_B1 = tx_ki_all_hla.DON_B1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_B2 = tx_ki_all_hla.DON_B2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_C1 = tx_ki_all_hla.DON_C1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_C2 = tx_ki_all_hla.DON_C2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DR1 = tx_ki_all_hla.DON_DR1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DR2 = tx_ki_all_hla.DON_DR2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DQ1 = tx_ki_all_hla.DON_DQ1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DQ2 = tx_ki_all_hla.DON_DQ2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DP1 = tx_ki_all_hla.DON_DP1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DP2 = tx_ki_all_hla.DON_DP2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DR51 = tx_ki_all_hla.DON_DR51.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DR52 = tx_ki_all_hla.DON_DR52.str.split(': ',expand=True)[1]
# tx_ki_all_hla.DON_DR53 = tx_ki_all_hla.DON_DR53.str.split(': ',expand=True)[1]

# tx_ki_all_hla.REC_A1 = tx_ki_all_hla.REC_A1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_A2 = tx_ki_all_hla.REC_A2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_B1 = tx_ki_all_hla.REC_B1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_B2 = tx_ki_all_hla.REC_B2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_CW1 = tx_ki_all_hla.REC_CW1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_CW2 = tx_ki_all_hla.REC_CW2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DR1 = tx_ki_all_hla.REC_DR1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DR2 = tx_ki_all_hla.REC_DR2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DQW1 = tx_ki_all_hla.REC_DQW1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DQW2 = tx_ki_all_hla.REC_DQW2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DPW1 = tx_ki_all_hla.REC_DPW1.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DPW2 = tx_ki_all_hla.REC_DPW2.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DRW51 = tx_ki_all_hla.REC_DRW51.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DRW52 = tx_ki_all_hla.REC_DRW52.str.split(': ',expand=True)[1]
# tx_ki_all_hla.REC_DRW53 = tx_ki_all_hla.REC_DRW53.str.split(': ',expand=True)[1]


# skip cases where no typing at A, B, or DR
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_A1'] != "Unknown"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_B1'] != "Unknown"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_DR1'] != "Unknown"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_A1'] != "Unknown"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_B1'] != "Unknown"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_DR1'] != "Unknown"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_A1'] != "Not Tested"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_B1'] != "Not Tested"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_DR1'] != "Not Tested"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_A1'] != "Not Tested"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_B1'] != "Not Tested"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_DR1'] != "Not Tested"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_A1'] != "0"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_B1'] != "0"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['DON_DR1'] != "0"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_A1'] != "0"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_B1'] != "0"]
tx_df_all_hla = tx_df_all_hla[tx_df_all_hla['REC_DR1'] != "0"]

print ("Minimum Typing of A, B, DR antigens: " + str(len(tx_df_all_hla)))

# missing DQ and DP - just skip the locus
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ1 == 'Unknown', 'DON_DQ1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == 'Unknown', 'DON_DP1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ1 == '0', 'DON_DQ1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == '0', 'DON_DP1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == 'Unknown', 'DON_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == 'Unknown', 'DON_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == '0', 'DON_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == '0', 'DON_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == 'Not Tested', 'DON_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == 'Not Tested', 'DON_DPA1'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_DQW1 == 'Unknown', 'REC_DQW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == 'Unknown', 'REC_DPW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW1 == '0', 'REC_DQW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == '0', 'REC_DPW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == 'Unknown', 'REC_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == 'Unknown', 'REC_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == '0', 'REC_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == '0', 'REC_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == 'Not Tested', 'REC_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == 'Not Tested', 'REC_DPA1'] = ""

# manage homozygosity

tx_df_all_hla.loc[tx_df_all_hla.DON_A2 == 'No second antigen detected', 'DON_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_B2 == 'No second antigen detected', 'DON_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DR2 == 'No second antigen detected', 'DON_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA2 == 'No second antigen detected', 'DON_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ2 == 'No second antigen detected', 'DON_DQ2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == 'No second antigen detected', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == 'No second antigen detected', 'DON_DP2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_A2 == 'Not Tested', 'DON_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_B2 == 'Not Tested', 'DON_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DR2 == 'Not Tested', 'DON_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA2 == 'Not Tested', 'DON_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ2 == 'Not Tested', 'DON_DQ2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == 'Not Tested', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == 'Not Tested', 'DON_DP2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_A2 == 'Unknown', 'DON_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_B2 == 'Unknown', 'DON_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DR2 == 'Unknown', 'DON_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA2 == 'Unknown', 'DON_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ2 == 'Unknown', 'DON_DQ2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == 'Unknown', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == 'Unknown', 'DON_DP2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_A2 == '0', 'DON_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_B2 == '0', 'DON_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DR2 == '0', 'DON_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA2 == '0', 'DON_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ2 == '0', 'DON_DQ2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == '0', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == '0', 'DON_DP2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_A2 == 'No second antigen detected', 'REC_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_B2 == 'No second antigen detected', 'REC_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DR2 == 'No second antigen detected', 'REC_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA2 == 'No second antigen detected', 'REC_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW2 == 'No second antigen detected', 'REC_DQW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == 'No second antigen detected', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == 'No second antigen detected', 'REC_DPW2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_A2 == 'Not Tested', 'REC_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_B2 == 'Not Tested', 'REC_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DR2 == 'Not Tested', 'REC_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA2 == 'Not Tested', 'REC_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW2 == 'Not Tested', 'REC_DQW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == 'Not Tested', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == 'Not Tested', 'REC_DPW2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_A2 == 'Unknown', 'REC_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_B2 == 'Unknown', 'REC_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DR2 == 'Unknown', 'REC_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA2 == 'Unknown', 'REC_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW2 == 'Unknown', 'REC_DQW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == 'Unknown', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == 'Unknown', 'REC_DPW2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_A2 == '0', 'REC_A2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_B2 == '0', 'REC_B2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DR2 == '0', 'REC_DR2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA2 == '0', 'REC_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW2 == '0', 'REC_DQW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == '0', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == '0', 'REC_DPW2'] = ""

# manage No 2nd antigen in first column - Some will result in two blanks but this is handled later

tx_df_all_hla.loc[tx_df_all_hla.DON_A1 == 'No second antigen detected', 'DON_A1'] = tx_df_all_hla.DON_A2
tx_df_all_hla.loc[tx_df_all_hla.DON_B1 == 'No second antigen detected', 'DON_B1'] = tx_df_all_hla.DON_B2
tx_df_all_hla.loc[tx_df_all_hla.DON_DR1 == 'No second antigen detected', 'DON_DR1'] = tx_df_all_hla.DON_DR2
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == 'No second antigen detected', 'DON_DQA1'] = tx_df_all_hla.DON_DQA2
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ1 == 'No second antigen detected', 'DON_DQ1'] = tx_df_all_hla.DON_DQ2
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == 'No second antigen detected', 'DON_DPA1'] = tx_df_all_hla.DON_DPA2
tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == 'No second antigen detected', 'DON_DP1'] = tx_df_all_hla.DON_DP2

tx_df_all_hla.loc[tx_df_all_hla.REC_A1 == 'No second antigen detected', 'REC_A1'] = tx_df_all_hla.REC_A2
tx_df_all_hla.loc[tx_df_all_hla.REC_B1 == 'No second antigen detected', 'REC_B1'] = tx_df_all_hla.REC_B2
tx_df_all_hla.loc[tx_df_all_hla.REC_DR1 == 'No second antigen detected', 'REC_DR1'] = tx_df_all_hla.REC_DR2
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == 'No second antigen detected', 'REC_DQA1'] = tx_df_all_hla.REC_DQA2
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW1 == 'No second antigen detected', 'REC_DQW1'] = tx_df_all_hla.REC_DQW2
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == 'No second antigen detected', 'REC_DPA1'] = tx_df_all_hla.REC_DPA2
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == 'No second antigen detected', 'REC_DPW1'] = tx_df_all_hla.REC_DPW2


# manage 4-digit typing without colons (DP has colons for 4-digit alelels)
 
# # https://stackoverflow.com/questions/47293729/replacing-part-of-a-string-in-pandas-column-with-and-condition
# regex=true for older Pandas version - 
tx_df_all_hla['DON_A1'] = tx_df_all_hla['DON_A1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_A2'] = tx_df_all_hla['DON_A2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_B1'] = tx_df_all_hla['DON_B1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_B2'] = tx_df_all_hla['DON_B2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_C1'] = tx_df_all_hla['DON_C1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_C2'] = tx_df_all_hla['DON_C2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_DR1'] = tx_df_all_hla['DON_DR1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_DR2'] = tx_df_all_hla['DON_DR2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_DQ1'] = tx_df_all_hla['DON_DQ1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['DON_DQ2'] = tx_df_all_hla['DON_DQ2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)


tx_df_all_hla['REC_A1'] = tx_df_all_hla['REC_A1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_A2'] = tx_df_all_hla['REC_A2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_B1'] = tx_df_all_hla['REC_B1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_B2'] = tx_df_all_hla['REC_B2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_CW1'] = tx_df_all_hla['REC_CW1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_CW2'] = tx_df_all_hla['REC_CW2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_DR1'] = tx_df_all_hla['REC_DR1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_DR2'] = tx_df_all_hla['REC_DR2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_DQW1'] = tx_df_all_hla['REC_DQW1'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)
tx_df_all_hla['REC_DQW2'] = tx_df_all_hla['REC_DQW2'].str.replace(r'^(\d\d)(\d\d)$',r'\1:\2', regex=True)


# manage C serology - No second antigen can mean low expression alleles - Imputation algorithm should do this!


# manage not tested and no second antigen for C and DQ and DQA1

tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == 'Not Tested', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == 'Not Tested', 'REC_CW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C2 == 'Not Tested', 'DON_C2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW2 == 'Not Tested', 'REC_CW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ1 == 'Not Tested', 'DON_DQ1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW1 == 'Not Tested', 'REC_DQW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ2 == 'Not Tested', 'DON_DQ2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW2 == 'Not Tested', 'REC_DQW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == 'Not Tested', 'DON_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == 'Not Tested', 'REC_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA2 == 'Not Tested', 'DON_DQA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA2 == 'Not Tested', 'REC_DQA2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == 'No antigen detected', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == 'No antigen detected', 'REC_CW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ1 == 'No antigen detected', 'DON_DQ1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW1 == 'No antigen detected', 'REC_DQW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == 'No antigen detected', 'DON_DQA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == 'No antigen detected', 'REC_DQA1'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == 'Unknown', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == 'No antigen detected', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == '0', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == 'Unknown', 'REC_CW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == 'No antigen detected', 'REC_CW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == '0', 'REC_CW1'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_C2 == 'Unknown', 'DON_C2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C2 == 'No second antigen detected', 'DON_C2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C2 == '0', 'DON_C2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW2 == 'Unknown', 'REC_CW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW2 == 'No second antigen detected', 'REC_CW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW2 == '0', 'REC_CW2'] = ""


tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == 'No second antigen detected', 'DON_C1'] = tx_df_all_hla.DON_C2
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == 'No second antigen detected', 'REC_CW1'] = tx_df_all_hla.REC_CW2
tx_df_all_hla.loc[tx_df_all_hla.DON_DQ1 == 'No second antigen detected', 'DON_DQ1'] = tx_df_all_hla.DON_DQ2
tx_df_all_hla.loc[tx_df_all_hla.REC_DQW1 == 'No second antigen detected', 'REC_DQW1'] = tx_df_all_hla.REC_DQW2
tx_df_all_hla.loc[tx_df_all_hla.DON_DQA1 == 'No second antigen detected', 'DON_DQA1'] = tx_df_all_hla.DON_DQA2
tx_df_all_hla.loc[tx_df_all_hla.REC_DQA1 == 'No second antigen detected', 'REC_DQA1'] = tx_df_all_hla.REC_DQA2

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

tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == '11', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C2 == '11', 'DON_C2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == '11', 'REC_CW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW2 == '11', 'REC_CW2'] = ""

# C13 - Invalid category - deleted allele

tx_df_all_hla.loc[tx_df_all_hla.DON_C1 == '13', 'DON_C1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_C2 == '13', 'DON_C2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW1 == '13', 'REC_CW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_CW2 == '13', 'REC_CW2'] = ""

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
tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == 'Not Tested', 'DON_DP1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == 'Not Tested', 'DON_DP2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == 'Not Tested', 'DON_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == 'Not Tested', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == 'Not Tested', 'REC_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == 'Not Tested', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == 'Not Tested', 'REC_DPW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == 'Not Tested', 'REC_DPW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == 'No second antigen detected', 'DON_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == 'No second antigen detected', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == 'No second antigen detected', 'DON_DP1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == 'No second antigen detected', 'DON_DP2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == 'No second antigen detected', 'REC_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == 'No second antigen detected', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == 'No second antigen detected', 'REC_DPW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == 'No second antigen detected', 'REC_DPW2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA1 == 'No antigen detected', 'DON_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DPA2 == 'No antigen detected', 'DON_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == 'No antigen detected', 'DON_DP1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == 'No antigen detected', 'DON_DP2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA1 == 'No antigen detected', 'REC_DPA1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPA2 == 'No antigen detected', 'REC_DPA2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == 'No antigen detected', 'REC_DPW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == 'No antigen detected', 'REC_DPW2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_DP1 == '1 (Inactive)', 'DON_DP1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DP2 == '1 (Inactive)', 'DON_DP2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW1 == '1 (Inactive)', 'REC_DPW1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DPW2 == '1 (Inactive)', 'REC_DPW2'] = ""


# Negative DRB3/4/5 typing - convert to blank and create DRB345 locus genotypes in next step
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB3_1 == 'N-Negative', 'REC_DRB3_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB3_2 == 'N-Negative', 'REC_DRB3_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB4_1 == 'N-Negative', 'REC_DRB4_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB4_2 == 'N-Negative', 'REC_DRB4_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB5_1 == 'N-Negative', 'REC_DRB5_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB5_2 == 'N-Negative', 'REC_DRB5_2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_DRB3_1 == 'N-Negative', 'DON_DRB3_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB3_2 == 'N-Negative', 'DON_DRB3_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB4_1 == 'N-Negative', 'DON_DRB4_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB4_2 == 'N-Negative', 'DON_DRB4_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB5_1 == 'N-Negative', 'DON_DRB5_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB5_2 == 'N-Negative', 'DON_DRB5_2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_DRB3_1 == '0', 'DON_DRB3_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB3_2 == '0', 'DON_DRB3_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB4_1 == '0', 'DON_DRB4_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB4_2 == '0', 'DON_DRB4_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB5_1 == '0', 'DON_DRB5_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB5_2 == '0', 'DON_DRB5_2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_DRB3_1 == '0', 'REC_DRB3_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB3_2 == '0', 'REC_DRB3_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB4_1 == '0', 'REC_DRB4_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB4_2 == '0', 'REC_DRB4_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB5_1 == '0', 'REC_DRB5_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB5_2 == '0', 'REC_DRB5_2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.DON_DRB3_1 == 'Not Tested', 'DON_DRB3_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB3_2 == 'Not Tested', 'DON_DRB3_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB4_1 == 'Not Tested', 'DON_DRB4_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB4_2 == 'Not Tested', 'DON_DRB4_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB5_1 == 'Not Tested', 'DON_DRB5_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.DON_DRB5_2 == 'Not Tested', 'DON_DRB5_2'] = ""

tx_df_all_hla.loc[tx_df_all_hla.REC_DRB3_1 == 'Not Tested', 'REC_DRB3_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB3_2 == 'Not Tested', 'REC_DRB3_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB4_1 == 'Not Tested', 'REC_DRB4_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB4_2 == 'Not Tested', 'REC_DRB4_2'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB5_1 == 'Not Tested', 'REC_DRB5_1'] = ""
tx_df_all_hla.loc[tx_df_all_hla.REC_DRB5_2 == 'Not Tested', 'REC_DRB5_2'] = ""


tx_df_all_hla.loc[(tx_df_all_hla.DON_DR51 == '95') & (tx_df_all_hla.DON_DRB5_1 == ''), 'DON_DRB5_1'] = "51"
tx_df_all_hla.loc[(tx_df_all_hla.DON_DR52 == '95') & (tx_df_all_hla.DON_DRB3_1 == ''), 'DON_DRB3_1'] = "52"
tx_df_all_hla.loc[(tx_df_all_hla.DON_DR53 == '95') & (tx_df_all_hla.DON_DRB4_1 == ''), 'DON_DRB4_1'] = "53"

tx_df_all_hla.loc[(tx_df_all_hla.REC_DRW51 == '95') & (tx_df_all_hla.REC_DRB5_1 == ''), 'REC_DRB5_1'] = "51"
tx_df_all_hla.loc[(tx_df_all_hla.REC_DRW52 == '95') & (tx_df_all_hla.REC_DRB3_1 == ''), 'REC_DRB3_1'] = "52"
tx_df_all_hla.loc[(tx_df_all_hla.REC_DRW53 == '95') & (tx_df_all_hla.REC_DRB4_1 == ''), 'REC_DRB4_1'] = "53"

# view Pandas data types
print (tx_df_all_hla.info())


# output merged HLA dataset

hla_filename = f"{donor_type}_{tx_type}_hla_9loc.csv"
tx_df_all_hla.to_csv(hla_filename, header=True, index=False)
print(f"Final {donor_type} {tx_type} HLA data saved to: {hla_filename}")

