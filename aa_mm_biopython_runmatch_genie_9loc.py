#!/usr/bin/env python

from pickle import FALSE
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import pandas as pd
import re
import gzip
import numpy as np
import os
import glob
from tqdm import tqdm
import hlagenie
genie = hlagenie.init("3510", ungap = True)
from aa_matching_msf_genie import *
aa_mm = AAMatch(dbversion=3420)

#loci = ["A","C","B","DRB1","DRB345","DQA1","DQB1","DPA1","DPB1"]

loci = ["A","C","B","DRB1","DQA1","DQB1","DPA1","DPB1"]


# set random seed
random.seed("20180413") # had to change seed to avoid deleted C*13:01

# flag to create runMatchMC files
generateRunMatchMC = 1



loci_value = 9
loci_range_value = 10

# SFVT file loading

# load variant type information file
SFVT_positions = {} # list of positions for each SFVT
SFVT_list = {} # list of SFVTs
seqf_name_ID = {} # sequence feature ID for each sequence feature name - to join to feature_2_name
#seqf_filename = "/Users/gracewager/dev/kidney-outcomes-sfvt/variant_type.bcp"
seqf_filename = "variant_type.bcp" 
seqf_file = open (seqf_filename,'r')
for line in seqf_file:
	# print (line)

	# skip whole alleles - seqf ID = 1 
	(iid, seqf_id, sfvt_name, *rest) = line.split('\t')
	if (not sfvt_name.startswith("Hsa")):
		continue

	(iid, seqf_id, sfvt_name, aa_string, SAP_list, ref_compare) = line.split('\t')
    # Hsa_HLA-A_SF220_VT474 becomes seqf_num "220"
	(hsa,hla_locus,seqf_num,vt_num) = sfvt_name.split("_")

	# skip complete protein
	if (hla_locus == ""): 
		continue
	
	# skip Unknown Motif Sequence
	if (vt_num == "Type Unknown"):
		continue

	seqf_num = re.sub("[A-Z]","",seqf_num) # remove "SF"
	(hla,locus) = hla_locus.split('-')  # determine locus

    # print (sfvt_name + " " + hla_locus + " " + seqf_num + " " + vt_num + "\n")

    # skip complete protein
	if (seqf_num == 1):
		continue

    # skip HLA-C SF 186 - missing data
	if ((locus == "C") & (seqf_num == "186")):
		continue

	# skip DQA1, DRA, DRB3, DRB4, DRB5, DPA1, DPB1 loci
	#if ((locus == "DQA1") | (locus == "DRA") | (locus == "DRB3") | (locus == "DRB4") | (locus == "DRB5") | (locus == "DPA1") | (locus == "DPB1")):
	#if ((locus == "DQA1") | (locus == "DRA") | (locus == "DPA1")):
		#continue
	if (locus == "DRA") | (locus == "DRBX"):
		continue

	# SAP.pm doesn't currently work outside of ARS
	# check that all polymorphic positions are in ARS
	SAPs = SAP_list.split("_")
	is_outside_ARS = 0
	positions = []
	for SAP in SAPs:
		position = SAP
		position = re.sub("[A-Za-z]","",position) # remove "ins" and "del"

		if ((locus == "A") | (locus == "C") | (locus == "B")):
			if ((int(position) < 1) | (int(position) > 183)):
				is_outside_ARS = 1
		else: # class II
			if ((int(position) < 5) | (int(position) > 95)):
				is_outside_ARS = 1
		# print (seqf_id + " " + position + "\n")
		positions.append(position)

	# skip SFVTs outside of ARS
	if (is_outside_ARS == 1):
		# print ("Seqf_ID outside of ARS: " + seqf_id)
		continue

	# skip single AAs because covered by SAP
	if (len(SAPs) == 1):
		continue 

	position_list = ','.join(positions)

	# 233 becomes "SFVT_A_233"
	seqf_name = "SFVT_" + locus + "_" + seqf_num

	seqf_name_ID[seqf_id] = seqf_name

	

A1_SFVT = []
A2_SFVT = []
C1_SFVT = []
C2_SFVT = []
B1_SFVT = []
B2_SFVT = []
DRB345_1_SFVT = []
DRB345_2_SFVT = []
DRB1_1_SFVT = []
DRB1_2_SFVT = []
DQA1_1_SFVT = []
DQA1_2_SFVT = []
DQB1_1_SFVT = []
DQB1_2_SFVT = []
DPA1_1_SFVT = []
DPA1_2_SFVT = []
DPB1_1_SFVT = []
DPB1_2_SFVT = []
for loc in loci:
    for id in range (1,1000):
        seqf_id = "SFVT_" + loc + "_" + str(id)
        if (seqf_id not in SFVT_list):
            continue
        if (loc == "A"):
            A1_SFVT_id = "A1_" + seqf_id
            A2_SFVT_id = "A2_" + seqf_id
            A1_SFVT.append(A1_SFVT_id)
            A2_SFVT.append(A2_SFVT_id)           
        elif (loc == "C"):
            C1_SFVT_id = "C1_" + seqf_id
            C2_SFVT_id = "C2_" + seqf_id
            C1_SFVT.append(C1_SFVT_id)
            C2_SFVT.append(C2_SFVT_id)
        elif (loc == "B"):
            B1_SFVT_id = "B1_" + seqf_id
            B2_SFVT_id = "B2_" + seqf_id
            B1_SFVT.append(B1_SFVT_id)
            B2_SFVT.append(B2_SFVT_id)
        elif (loc == "DRB3"):
            DRB345_1_SFVT_id = "DRB345_1_" + seqf_id
            DRB345_2_SFVT_id = "DRB345_2_" + seqf_id
            DRB345_1_SFVT.append(DRB345_1_SFVT_id)
            DRB345_2_SFVT.append(DRB345_2_SFVT_id)
        elif (loc == "DRB4"):
            DRB345_1_SFVT_id = "DRB345_1_" + seqf_id
            DRB345_2_SFVT_id = "DRB345_2_" + seqf_id
            DRB345_1_SFVT.append(DRB345_1_SFVT_id)
            DRB345_2_SFVT.append(DRB345_2_SFVT_id)
        elif (loc == "DRB5"):
            DRB345_1_SFVT_id = "DRB345_1_" + seqf_id
            DRB345_2_SFVT_id = "DRB345_2_" + seqf_id
            DRB345_1_SFVT.append(DRB345_1_SFVT_id)
            DRB345_2_SFVT.append(DRB345_2_SFVT_id)
        elif (loc == "DRB1"):
            DRB1_1_SFVT_id = "DRB1_1_" + seqf_id
            DRB1_2_SFVT_id = "DRB1_2_" + seqf_id
            DRB1_1_SFVT.append(DRB1_1_SFVT_id)
            DRB1_2_SFVT.append(DRB1_2_SFVT_id)
        elif (loc == "DQA1"):
            DQA1_1_SFVT_id = "DQA1_1_" + seqf_id
            DQA1_2_SFVT_id = "DQA1_2_" + seqf_id
            DQA1_1_SFVT.append(DQA1_1_SFVT_id)
            DQA1_2_SFVT.append(DQA1_2_SFVT_id)
        elif (loc == "DQB1"):
            DQB1_1_SFVT_id = "DQB1_1_" + seqf_id
            DQB1_2_SFVT_id = "DQB1_2_" + seqf_id
            DQB1_1_SFVT.append(DQB1_1_SFVT_id)
            DQB1_2_SFVT.append(DQB1_2_SFVT_id)
        elif (loc == "DPA1"):
            DPA1_1_SFVT_id = "DPA1_1_" + seqf_id
            DPA1_2_SFVT_id = "DPA1_2_" + seqf_id
            DPA1_1_SFVT.append(DPA1_1_SFVT_id)
            DPA1_2_SFVT.append(DPA1_2_SFVT_id)
        elif (loc == "DPB1"):
            DPB1_1_SFVT_id = "DPB1_1_" + seqf_id
            DPB1_2_SFVT_id = "DPB1_2_" + seqf_id
            DPB1_1_SFVT.append(DPB1_1_SFVT_id)
            DPB1_2_SFVT.append(DPB1_2_SFVT_id)			
        else: # data from other loci not currently used
            continue


A1_SFVT_header = ",".join(A1_SFVT)
A2_SFVT_header = ",".join(A2_SFVT)
C1_SFVT_header = ",".join(C1_SFVT)
C2_SFVT_header = ",".join(C2_SFVT)
B1_SFVT_header = ",".join(B1_SFVT)
B2_SFVT_header = ",".join(B2_SFVT)
DRB345_1_SFVT_header =  ",".join(DRB345_1_SFVT)
DRB345_2_SFVT_header =  ",".join(DRB345_2_SFVT)
DRB1_1_SFVT_header =  ",".join(DRB1_1_SFVT)
DRB1_2_SFVT_header =  ",".join(DRB1_2_SFVT)
DQA1_1_SFVT_header =  ",".join(DQA1_1_SFVT)
DQA1_2_SFVT_header =  ",".join(DQA1_2_SFVT)
DQB1_1_SFVT_header =  ",".join(DQB1_1_SFVT)
DQB1_2_SFVT_header =  ",".join(DQB1_2_SFVT)
DPA1_1_SFVT_header =  ",".join(DPA1_1_SFVT)
DPA1_2_SFVT_header =  ",".join(DPA1_2_SFVT)
DPB1_1_SFVT_header =  ",".join(DPB1_1_SFVT)
DPB1_2_SFVT_header =  ",".join(DPB1_2_SFVT)


# header_AA_string = ",".join([A1_AA_string,A2_AA_string,C1_AA_string,C2_AA_string,B1_AA_string,B2_AA_string,DRB1_1_AA_string,DRB1_2_AA_string,DQB1_1_AA_string,DQB1_2_AA_string])

header_SFVT_string = ','.join([A1_SFVT_header,A2_SFVT_header,C1_SFVT_header,C2_SFVT_header,B1_SFVT_header,B2_SFVT_header,DRB1_1_SFVT_header,DRB1_2_SFVT_header,DQB1_1_SFVT_header,DQB1_2_SFVT_header])

# SFVT mappings
sfvt_filename = "permissive_SFVT.txt"
sfvt_file = open (sfvt_filename,'w')
# sfvt_file.write(header_AA_string + "," + header_SFVT_string + "\n")
sfvt_file.write(header_SFVT_string + "\n")


# mismatch_list.txt from HSCT permissive study had up to 3 mismatched loci - solid organ can all loci mismatched

crids_MM = {} # crids from MM file
surv1yr= {} # survival at 1 year for CRID
agvhd= {} # acute GVHD for CRID
SAP_directional= {} # directional differences in amino acids
SAP_pair= {} # nondirectional differences in amino acids
SAP_position= {} # positions where alleles differ
SAP_directional_ID= {} # directional differences in amino acids
SAP_pair_ID= {} # nondirectional differences in amino acids
SAP_position_ID= {} # positions where alleles differ
SAP_directional_list= {} # directional differences in amino acids
SAP_pair_list= {} # nondirectional differences in amino acids
SAP_position_list= {} # positions where alleles differ
SFVT_directional= {} # directional differences in amino acids
SFVT_pair= {} # nondirectional differences in amino acids
SFVT_position= {} # positions where alleles differ
SFVT_directional_ID= {} # directional differences in amino acids
SFVT_pair_ID= {} # nondirectional differences in amino acids
SFVT_position_ID= {} # positions where alleles differ
SFVT_directional_list= {} # directional differences in amino acids
SFVT_pair_list= {} # nondirectional differences in amino acids
SFVT_position_list= {} # positions where alleles diffe


# TODO - from HLA transplant pair dataset, create array of mismatched allele pairs
# TODO - SAP_position - list mismatched positions per pair 
# TODO - SAP_directional - list actual AAs in donor and recip listing donor AA list A_A_
# TODO - SAP pair - list actual AAS in donor and recip in sorted order




# testing the functions

allele1 = "A*02:01"
allele2 = "A*01:01"


print (allele1)
print(type(allele1))
print(type(44))
AA = genie.getAA(allele1,6)
print ("AA at Position 6: " + AA)

#AA_string = genie.seqs(allele1)
#print ("AAs for mature protein: " + str(AA_string))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"A*02:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))

print("AAs for mature protein: ", len(AA_string_mature['A*02:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"B*07:02"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['B*07:02']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"C*01:02"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['C*01:02']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DRB1*01:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DRB1*01:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DRB3*01:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DRB3*01:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DRB4*01:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DRB4*01:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DRB5*01:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DRB5*01:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DQA1*01:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DQA1*01:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DQB1*05:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DQB1*05:01']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DPA1*01:03"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DPA1*01:03']))

AA_string_mature = {key: genie.seqs[key] for key in genie.seqs.keys()
       & {"DPB1*01:01"}}
print ("AAs for mature protein: " + str(AA_string_mature))
print("AAs for mature protein: ", len(AA_string_mature['DPB1*01:01']))

AA_string_full = {key: genie.full_seqs[key] for key in genie.full_seqs.keys()
       & {"A*02:01"}}
print ("AAs for full protein: " + str(AA_string_full))

#peptide_string = genie.getPeptide(allele1)
#print ("Peptide String: " + peptide_string)


ARD_string = genie.getARD(allele1)
print ("ARD: " + ARD_string)

epitope_string = genie.getEpitope(allele1,[1,3,5,8,10])
print ("Epitope for AA Position List [1,3,5,8,10]: " + epitope_string)

AA = genie.getAA(allele1,44)
print ("A*02:01 position 44: " + AA)
AA = genie.getAA(allele2,44)
print ("A*01:01 position 44: " + AA)

mismatched = genie.isPositionMismatched(allele1,allele2,44)
print ("Are A*02:01 and A*01:01 mismatched at position 44?: " + str(mismatched))


AA = genie.getAA(allele1,45)
print ("A*02:01 position 45: " + AA)
AA = genie.getAA(allele2,45)
print ("A*01:01 position 45: " + AA)

mismatched = genie.isPositionMismatched(allele1,allele2,45)
print ("Are A*02:01 and A*01:01 mismatched at position 45?: " + str(mismatched))

mm_count = genie.countAAMismatchesAllele(allele1,allele1,allele2,allele2,44)
print ("# of Mismatches between A*02:01 homozygous donor and A*01:01 homozygous recip at position 44: " + str(mm_count))

AA1_donor = "Y"
AA2_donor = "Y"
AA1_recip = "D"
AA2_recip = "D"
#mm_count = genie.countAAMismatchesAllele(AA1_donor,AA2_donor,AA1_recip,AA2_recip)
#print ("# of Mismatches between YY AAs in donor and DD in recip: " + str(mm_count))
#AA1_recip = "Y"
#mm_count = genie.countAAMismatchesAllele(AA1_donor,AA2_donor,AA1_recip,AA2_recip)
#print ("# of Mismatches between YY AAs in donor and YD in recip: " + str(mm_count))

allele1_donor = "A*02:01"
allele2_donor = "A*01:01"
allele1_recip = "A*03:01"
allele2_recip = "A*66:01"

# loop through all the SFVT_IDs
for id in range(1,1000):
	locus = "A"
	seqf_name= "SFVT_" + locus + "_"  + str(id)
	if (seqf_name not in SFVT_list):
		continue	
	
	position_list = SFVT_positions[seqf_name]
	positions = [int(x) for x in position_list.split(",")]

	SFVT_MM_count = aa_mm.count_AA_Mismatches_SFVT(allele1_donor,allele2_donor,allele1_recip,allele2_recip,positions)

	# print ("Sequence Feature Name: " + seqf_name + " Positions: " + position_list + " MM Count: " + str(SFVT_MM_count))

# load SRTR HapLogic imputation output file and load probabilities for all subjects
#all of these are 5 loci#
#impute_outfilename = pathloc + "/srtr_impute_10000.stage"
#impute_outfile = open(impute_outfilename, "rt")
#impute_outfilename = pathloc + "/srtr_impute.stage.gz"
#impute_outfile = gzip.open(impute_outfilename, "rt")

#join all ethnic data into one table
#haplo = os.path.join("./kamoun_impute/", "impute_srtr.*.csv")
#haplo = glob.glob(haplo)
#haplo = pd.concat(map(pd.read_csv, haplo), ignore_index=FALSE)
#haplo_matrix_filename = pathloc + "/impute_all_ethnicitys.csv"
#haplo.to_csv(haplo_matrix_filename, index= FALSE)


# compute cumulative genotype frequency totals per subject
happair_id_total = {} # cumulative genotype frequency total
happair_probs = {} # HLA probability distribution
happair_hla = {} # HLA haplotype pair distribution

# load SRTR HapLogic imputation output file and load probabilities for all subjects
pops = ['AFA','ASN','CAU','HIS','NAM','MLT']

for pop in pops:
	srtr_impute_filename = "impute.srtr." + pop + ".csv.gz"
	with gzip.open(srtr_impute_filename, 'rt') as impute_outfile:
		for line in impute_outfile:
			#print(line)
			(subject_id,cohort_id,hap1,hap2,freq) = line.split(',')
			if (cohort_id == "cohort_id"): # skip header row
				continue

			happair_freq = float(freq)

			if subject_id not in happair_id_total:
				happair_id_total[subject_id] = 0
			happair_id_total[subject_id] += float(happair_freq)


	impute_outfile.close()


# Single genotype frequnecy already calculated for haplo pairs in 9loc impute data 
for pop in pops:
	srtr_impute_filename = "impute.srtr." + pop + ".csv.gz"
	with gzip.open(srtr_impute_filename, 'rt') as impute_outfile:
		for line in impute_outfile:
			(subject_id,cohort_id,hap1,hap2,freq) = line.split(',')
			if (cohort_id == "cohort_id"): # skip header row
				continue
			#print(subject_id)    
			#print(line)
			happair_freq = 0
			happair_freq = float(freq)

			happair = hap1 + "+" + hap2

			if subject_id not in happair_probs:
				happair_probs[subject_id] = []
			freq_tots = happair_freq / happair_id_total[subject_id]
			happair_probs[subject_id].append(freq_tots)

			if subject_id not in happair_hla:
				happair_hla[subject_id]= []
			happair_hla[subject_id].append(happair)
	impute_outfile.close()


# load in outcomes data from TX_KI file (tab-delimited is best?)
tx_ki_filename = "TX_KI_decoded.txt"
tx_ki = pd.read_csv(tx_ki_filename, sep='\t', index_col=None, encoding='Latin-1', low_memory=False)

# subset TX_KI columns to those used in previous study
# TFL_ENDX became TFL_ENDTXFU in SRTR SAF since Keith's script
tx_ki = tx_ki[["TX_ID","PX_ID","REC_TX_DT","ORG_TY","PERS_ID","TFL_DEATH_DT","REC_AGE_AT_TX",
						"REC_TX_ORG_TY","CAN_GENDER","CAN_DGN","CAN_RACE","REC_COLD_ISCH_TM","CAN_AGE_DIAB",
						"CAN_DIAB","CAN_DIAB_TY","CAN_ABO","REC_PREV_KI","REC_PREV_KP","REC_HGT_CM",
						"REC_WGT_KG","REC_A_MM_EQUIV_CUR","REC_B_MM_EQUIV_CUR","REC_DR_MM_EQUIV_CUR",
						"DON_TY","DON_AGE","DON_CAD_DON_COD","DON_GENDER","DON_RACE","DON_EXPAND_DON_KI",
						"DON_RT_KI_PUMP","DON_LF_KI_PUMP","TFL_ENDTXFU","DON_CREAT","REC_MED_COND","REC_DGN",
						"DON_KI_CREAT_PREOP","DON_HIST_HYPERTEN","DON_RELATIONSHIP_TY","TFL_LAFUDATE",
						"PERS_RETX","REC_PRIMARY_PAY","REC_SECONDARY_PAY","REC_TX_PROCEDURE_TY",
						"DON_HGT_CM","DON_WGT_KG","DON_ORG_SHARED","DON_NON_HR_BEAT","DON_HTN",
						"REC_PRETX_TXFUS","DON_HIST_DIAB","REC_CMV_IGG","REC_CMV_IGM","DON_CMV_IGG",
						"DON_ANTI_CMV","CAN_MALIG","TFL_GRAFT_DT","CAN_LAST_SRTR_PEAK_PRA","REC_DIAL_DT",
						"DON_RELATIONSHIP_TY","CAN_SOURCE","REC_CREAT",
						"REC_PREV_PREG","REC_FIRST_WEEK_DIAL","DONOR_ID"]]
# print (tx_ki)

# load SRTR transplant pairs with both subject IDs

impute_err_filename = "impute_failed_ids.csv"
impute_err_file = open(impute_err_filename, "w")
impute_err_file.write("Type,ID\n")

#  loop through 10 replicate random realizations for multiple imputation
multiple_imputation_replicates = 10
for rep in tqdm(range(1,multiple_imputation_replicates+1)):
	print(rep)
	# tx_ki_hla.csv
	# ORG_TY,PERS_ID,PX_ID,REC_TX_DT,REC_HISTO_TX_ID,DON_TY,DON_RACE,DON_RACE_SRTR,DON_ETHNICITY_SRTR,DON_A1,DON_A2,DON_B1,DON_B2,DON_DR1,DON_DR2,REC_AGE_IN_MONTHS_AT_TX,CAN_RACE,CAN_RACE_SRTR,CAN_ETHNICITY_SRTR,REC_TX_TY,REC_A1,REC_A2,REC_B1,REC_B2,REC_DR1,REC_DR2,DONOR_ID,DON_C1,DON_C2,DON_DQ1,DON_DQ2,DON_DP1,DON_DP2,DON_DR51,DON_DR52,DON_DR53,REC_CW1,REC_CW2,REC_DPW1,REC_DPW2,REC_DQW1,REC_DQW2,REC_DRW51,REC_DRW52,REC_DRW53
	tx_hla_filename = "tx_ki_hla_9loc.csv"
	tx_hla_file = open(tx_hla_filename, "r")
	runmatch_filename = "out.runmatchMCgenie." + str(rep) + ".txt.gz"
	runmatch_file = gzip.open(runmatch_filename, "wt")

	#happair_hla_filename = pathloc + "/happair_hla.csv"
	#happair_hla_file = open(happair_hla_filename, "r")
	#print(happair_hla)

	happair_selected_donor = {} # haplotype pair selected by weighted choice
	happair_selected_recip = {} # haplotype pair selected by weighted choice
	loci_pos = []
	for line in tx_hla_file:
		(ORG_TY,PERS_ID,PX_ID,REC_TX_DT,REC_HISTO_TX_ID,DON_TY,DON_RACE,DON_RACE_SRTR,DON_ETHNICITY_SRTR,\
   			DON_A1,DON_A2,DON_B1,DON_B2,DON_DR1,DON_DR2,\
			REC_AGE_IN_MONTHS_AT_TX,CAN_RACE,CAN_RACE_SRTR,CAN_ETHNICITY_SRTR,REC_TX_TY,\
			REC_A1,REC_A2,REC_B1,REC_B2,REC_DR1,REC_DR2,\
			DONOR_ID,DON_C1,DON_C2,DON_DQ1,DON_DQ2,DON_DP1,DON_DP2,\
			DON_DR51,DON_DR52,DON_DR53,\
			REC_CW1,REC_CW2,REC_DPW1,REC_DPW2,REC_DQW1,REC_DQW2,REC_DRW51,REC_DRW52,REC_DRW53,\
			DON_DRB3_1,DON_DRB3_2,DON_DRB4_1,DON_DRB4_2,DON_DRB5_1,DON_DRB5_2,DON_DQA1,DON_DQA2,DON_DPA1,DON_DPA2,\
			REC_DRB3_1,REC_DRB3_2,REC_DRB4_1,REC_DRB4_2,REC_DRB5_1,REC_DRB5_2,REC_DQA1,REC_DQA2,REC_DPA1,REC_DPA2) = line.split(',')
		if (ORG_TY == "ORG_TY"):
			continue

		# DONOR - DONOR_ID
		# RECIP - PERS_ID
		# TRANSPLANT PAIR - PX_ID

		SUBJECT_ID_RECIP = "R" + PX_ID
		#print(SUBJECT_ID_RECIP)
		SUBJECT_ID_DONOR = "D" + PX_ID

		# check for missing imputation output
		if (SUBJECT_ID_RECIP not in happair_hla):
			# print ("Missing Recip ID: " + PERS_ID)
			if (rep == 1):
				impute_err_file.write(SUBJECT_ID_RECIP + "\n")
			continue

		if (SUBJECT_ID_DONOR not in happair_hla):
			# print ("Missing Donor ID: " + DONOR_ID)
			if (rep == 1):
				impute_err_file.write(SUBJECT_ID_DONOR + "\n")
			continue
		
		# select haplotype pair based on probability distribution
		haplistR = happair_hla[SUBJECT_ID_RECIP]
		problistR = happair_probs[SUBJECT_ID_RECIP]

		happair_recip = aa_mm.weighted_choice(haplistR,problistR)
		#print ("Random Recip HapPair: " + PERS_ID + " " + happair_recip)

		haplistD = happair_hla[SUBJECT_ID_DONOR]
		problistD = happair_probs[SUBJECT_ID_DONOR]

		happair_donor = aa_mm.weighted_choice(haplistD,problistD)
		#print ("Random Donor HapPair: " + DONOR_ID + " " + happair_donor)
		(hap1_donor,hap2_donor) = happair_donor.split('+')
		(a1_donor,c1_donor,b1_donor,drb345_1_donor,drb1_1_donor,dqa1_1_donor,dqb1_1_donor,dpa1_1_donor,dpb1_1_donor) = hap1_donor.split('~')
		(a2_donor,c2_donor,b2_donor,drb345_2_donor,drb1_2_donor,dqa1_2_donor,dqb1_2_donor,dpa1_2_donor,dpb1_2_donor) = hap2_donor.split('~')

		(hap1_recip,hap2_recip) = happair_recip.split('+')
		(a1_recip,c1_recip,b1_recip,drb345_1_recip,drb1_1_recip,dqa1_1_recip,dqb1_1_recip,dpa1_1_recip,dpb1_1_recip) = hap1_recip.split('~')
		(a2_recip,c2_recip,b2_recip,drb345_2_recip,drb1_2_recip,dqa1_2_recip,dqb1_2_recip,dpa1_2_recip,dpb1_2_recip) = hap2_recip.split('~')

		# handle deleted / renamed alleles between IMGT/HLA 3.4 and 3.43

		if (a1_donor == "A*23:19Q"):
			a1_donor = "A*23:19N"
		if (a2_donor == "A*23:19Q"):
			a2_donor = "A*23:19N"
		if (a1_recip == "A*23:19Q"):
			a1_recip = "A*23:19N"
		if (a2_recip == "A*23:19Q"):
			a2_recip = "A*23:19N"

		if (b1_donor == "B*15:22"):
			b1_donor = "B*35:43"
		if (b2_donor == "B*15:22"):
			b2_donor = "B*35:43"
		if (b1_recip == "B*15:22"):
			b1_recip = "B*35:43"
		if (b2_recip == "B*15:22"):
			b2_recip = "B*35:43"

		if (c1_donor == "C*15:20"):
			c1_donor = "C*15:27"
		if (c2_donor == "C*15:20"):
			c2_donor = "C*15:27"
		if (c1_recip == "C*15:20"):
			c1_recip = "C*15:27"
		if (c2_recip == "C*15:20"):
			c2_recip = "C*15:27"

		if (c1_donor == "C*05:02"):
			c1_donor = "C*05:09"
		if (c2_donor == "C*05:02"):
			c2_donor = "C*05:09"
		if (c1_recip == "C*05:02"):
			c1_recip = "C*05:09"
		if (c2_recip == "C*05:02"):
			c2_recip = "C*05:09"

		if (c1_donor == "C*03:12"):
			c1_donor = "C*03:19"
		if (c2_donor == "C*03:12"):
			c2_donor = "C*03:19"
		if (c1_recip == "C*03:12"):
			c1_recip = "C*03:19"
		if (c2_recip == "C*03:12"):
			c2_recip = "C*03:19"

		if (c1_donor == "C*03:23"):
			c1_donor = "C*03:23N"
		if (c2_donor == "C*03:23"):
			c2_donor = "C*03:23N"
		if (c1_recip == "C*03:23"):
			c1_recip = "C*03:23N"
		if (c2_recip == "C*03:23"):
			c2_recip = "C*03:23N"

		if (drb1_1_donor == "DRB1*07:02"):
			drb1_1_donor = "DRB1*07:01"
		if (drb1_2_donor == "DRB1*07:02"):
			drb1_2_donor = "DRB1*07:01"
		if (drb1_1_recip == "DRB1*07:02"):
			drb1_1_recip = "DRB1*07:01"
		if (drb1_2_recip == "DRB1*07:02"):
			drb1_2_recip = "DRB1*07:01"

		# C*13:01 was deleted long ago - not sure changing to 12:02 is the right move - imputation algorithm issue
		if (c1_donor == "C*13:01"):
			c1_donor = "C*12:02"
		if (c2_donor == "C*13:01"):
			c2_donor = "C*12:02"
		if (c1_recip == "C*13:01"):
			c1_recip = "C*12:02"
		if (c2_recip == "C*13:01"):
			c2_recip = "C*12:02"


		if (dqa1_1_donor == "DQA1*01:07"):
			dqa1_1_donor = "DQA1*01:07Q"
		if (dqa1_2_donor == "DQA1*01:07"):
			dqa1_2_donor = "DQA1*01:07Q"
		if (dqa1_1_recip == "DQA1*01:07"):
			dqa1_1_recip = "DQA1*01:07Q"
		if (dqa1_2_recip == "DQA1*01:07"):
			dqa1_2_recip = "DQA1*01:07Q"


		# store chosen haplos based on PX_ID to merge into TX_KI later
		hap1_donor_alleles = [a1_donor,c1_donor,b1_donor,drb345_1_donor,drb1_1_donor,dqa1_1_donor,dqb1_1_donor,dpa1_1_donor,dpb1_1_donor]
		hap1_donor = '~'.join(hap1_donor_alleles)
		hap2_donor_alleles = [a2_donor,c2_donor,b2_donor,drb345_2_donor,drb1_2_donor,dqa1_2_donor,dqb1_2_donor,dpa1_2_donor,dpb1_2_donor]
		hap2_donor = '~'.join(hap2_donor_alleles)
		happair_donor = hap1_donor + "+" + hap2_donor
		hap1_recip_alleles = [a1_recip,c1_recip,b1_recip,drb345_1_recip,drb1_1_recip,dqa1_1_recip,dqb1_1_recip,dpa1_1_recip,dpb1_1_recip]
		hap1_recip = '~'.join(hap1_recip_alleles)
		hap2_recip_alleles = [a2_recip,c2_recip,b2_recip,drb345_2_recip,drb1_2_recip,dqa1_2_recip,dqb1_2_recip,dpa1_2_recip,dpb1_2_recip]
		hap2_recip = '~'.join(hap2_recip_alleles)
		happair_recip = hap1_recip + "+" + hap2_recip
		#print(happair_recip)

		happair_selected_recip[PX_ID] = happair_recip
		happair_selected_donor[PX_ID] = happair_donor
		# happair_selected_recip_df.append([PX_ID,happair_recip])
		# happair_selected_donor_df.append([PX_ID,happair_donor])

		# stop here if not writing runmatch files
		if (generateRunMatchMC == 0):
			continue

		# runmatchMC file


		#print ("Multiple imputation replicate: " + runmatch_filename)

		for loc in loci:
			#print(loc)

			#handeling of DRBX
			if (loc == "A"):
				allele1_donor = a1_donor
				allele2_donor = a2_donor
				allele1_recip = a1_recip
				allele2_recip = a2_recip
			if (loc == "C"):
				allele1_donor = c1_donor
				allele2_donor = c2_donor
				allele1_recip = c1_recip
				allele2_recip = c2_recip
			if (loc == "B"):
				allele1_donor = b1_donor
				allele2_donor = b2_donor
				allele1_recip = b1_recip
				allele2_recip = b2_recip
			if (loc == "DRB345"):
				allele1_donor = drb345_1_donor
				allele2_donor = drb345_2_donor
				allele1_recip = drb345_1_recip
				allele2_recip = drb345_2_recip
			if (loc == "DRB1"):
				allele1_donor = drb1_1_donor
				allele2_donor = drb1_2_donor
				allele1_recip = drb1_1_recip
				allele2_recip = drb1_2_recip
			if (loc == "DQA1"):
				allele1_donor = dqa1_1_donor
				allele2_donor = dqa1_2_donor
				allele1_recip = dqa1_1_recip
				allele2_recip = dqa1_2_recip
			if (loc == "DQB1"):
				allele1_donor = dqb1_1_donor
				allele2_donor = dqb1_2_donor
				allele1_recip = dqb1_1_recip
				allele2_recip = dqb1_2_recip
			if (loc == "DPA1"):
				allele1_donor = dpa1_1_donor
				allele2_donor = dpa1_2_donor
				allele1_recip = dpa1_1_recip
				allele2_recip = dpa1_2_recip
			if (loc == "DPB1"):
				allele1_donor = dpb1_1_donor
				allele2_donor = dpb1_2_donor
				allele1_recip = dpb1_1_recip
				allele2_recip = dpb1_2_recip

			full_start = aa_mm.full_start_pos[loc]
			full_end = aa_mm.full_end_pos[loc]

			# count mismatches for positions where there will be complete sequences
			for pos in range(full_start,full_end):
				# AA = getAAposition(HLA_seq,allele1_donor,pos)
				# print (allele1_donor + " - " + str(pos) + AA)
				mm_count = genie.countAAMismatchesAllele(allele1_donor,allele2_donor,allele1_recip,allele2_recip,pos)
				#print ("MM_count at position " + str(pos) + ": " + str(mm_count) + " - HLA - " + allele1_donor + " " + allele2_donor + " " + allele1_recip + " " + allele2_recip)

				mm_count_0 = 0
				mm_count_1 = 0
				mm_count_2 = 0
				if (mm_count == 0):
					mm_count_0 = 1
				if (mm_count == 1):
					mm_count_1 = 1
				if (mm_count == 2):
					mm_count_2 = 1


				#print (PX_ID + "|" + loc + "|" + str(pos) + "|" + str(mm_count_0) + "|" + str(mm_count_1) + "|" + str(mm_count_2))

				# create 5-locus out.runmatchMC files on the fly
				runmatch_file.write(PX_ID + "|" + loc + "|" + str(pos) + "|" + str(mm_count_0) + "|" + str(mm_count_1) + "|" + str(mm_count_2) + "\n")
	tx_hla_file.close()

	# DATA MATRIX FILE OUTPUT

	# create dataframes from happair dictionaries
	happair_selected_recip_df = pd.DataFrame(list(happair_selected_recip.items()),columns=["PX_ID","HAPPAIR_RECIP"])
	happair_selected_donor_df = pd.DataFrame(list(happair_selected_donor.items()),columns=["PX_ID","HAPPAIR_DONOR"])

	happair_selected = happair_selected_recip_df.merge(happair_selected_donor_df,left_on="PX_ID",right_on="PX_ID")
	# happair_selected = pd.concat(happair_selected_donor_df,happair_selected_recip_df,axis=1)
	'''
	with pd.option_context("display.max_rows",None,"display.max_columns",None):
		print (happair_selected)
	'''

	# create dataframes with additional HLA columns
	#print ("selected:", happair_selected_recip_df["HAPPAIR_RECIP"])

	# new data frame with split value columns for loci
	recip_haplos = happair_selected["HAPPAIR_RECIP"].str.split("+", n = 1, expand = True)
	donor_haplos = happair_selected["HAPPAIR_DONOR"].str.split("+", n = 1, expand = True)
	happair_selected["RECIP_HAP1"] = recip_haplos[0]
	happair_selected["RECIP_HAP2"] = recip_haplos[1]
	happair_selected["DONOR_HAP1"] = donor_haplos[0]
	happair_selected["DONOR_HAP2"] = donor_haplos[1]
	recip_alleles1 = happair_selected["RECIP_HAP1"].str.split("~", n = loci_value, expand = True)
	recip_alleles2 = happair_selected["RECIP_HAP2"].str.split("~", n = loci_value, expand = True)
	donor_alleles1 = happair_selected["DONOR_HAP1"].str.split("~", n = loci_value, expand = True)
	donor_alleles2 = happair_selected["DONOR_HAP2"].str.split("~", n = loci_value, expand = True)
	# happair_selected["RECIP_A1"] = recip_alleles1[0]
	# print (happair_selected)
	# print (recip_alleles1)

	# add allele columns
	for i in range(1,loci_range_value):
		loc = loci[i-1]
		recip_alleles1_column_name = "RECIP_" + loc + "_1"
		recip_alleles2_column_name = "RECIP_" + loc + "_2"
		donor_alleles1_column_name = "DONOR_" + loc + "_1"
		donor_alleles2_column_name = "DONOR_" + loc + "_2"
		recip_allele1 = recip_alleles1[i-1]
		recip_allele2 = recip_alleles2[i-1]
		donor_allele1 = donor_alleles1[i-1]
		donor_allele2 = donor_alleles2[i-1]		
		happair_selected[recip_alleles1_column_name] = recip_allele1
		happair_selected[recip_alleles2_column_name] = recip_allele2
		happair_selected[donor_alleles1_column_name] = donor_allele1
		happair_selected[donor_alleles2_column_name] = donor_allele2

	# add AA position columns
	# for i in range(1,5):
	#	loc = loci[i-1]
	# 	ard_start = ard_start_pos[loc]
	# 	ard_end = ard_end_pos[loc]
	# 	for pos in range(ard_start,ard_end):
	# 		recip_aa1_column_name = "RECIP_" + loc + "_1_" + str(pos)
	# 		recip_aa2_column_name = "RECIP_" + loc + "_2_" + str(pos)
	# 		donor_aa1_column_name = "DONOR_" + loc + "_1_" + str(pos)
	# 		donor_aa2_column_name = "DONOR_" + loc + "_2_" + str(pos)

			# initialize new columns
	# 		happair_selected[recip_aa1_column_name] = happair_selected[recip_alleles1_column_name]
	# 		happair_selected[recip_aa2_column_name] = happair_selected[recip_alleles2_column_name]
	# 		happair_selected[donor_aa1_column_name] = happair_selected[donor_alleles1_column_name]
	# 		happair_selected[donor_aa2_column_name] = happair_selected[donor_alleles2_column_name]

			# get AA position and count mismatches for each row
	# 		for index,row in happair_selected.head().iterrows():
	# 			happair_selected[index,recip_aa1_column_name] = getAAposition(HLA_seq,recip_allele1[index],pos)
	# 			happair_selected[index,recip_aa2_column_name] = getAAposition(HLA_seq,recip_allele2[index],pos)
	# 			happair_selected[index,donor_aa1_column_name] = getAAposition(HLA_seq,donor_allele1[index],pos)
	# 			happair_selected[index,donor_aa2_column_name] = getAAposition(HLA_seq,donor_allele2[index],pos)

	# add AA mismatch columns
	for i in range(1,loci_range_value):
		recip_allele1 = recip_alleles1[i-1]
		loc = recip_allele1[0].split('*')[0]
		print(loc)
		if loc in ['DRB3', 'DRB4', 'DRB5', 'DRBX']:
			loc = 'DRB345'
			full_start = aa_mm.full_start_pos['DRB3']
			full_end = aa_mm.full_end_pos['DRB3']
		else:
			full_start = aa_mm.full_start_pos[loc]
			full_end = aa_mm.full_end_pos[loc]		
		#print("recip_alleles1: "+str(recip_alleles1)+"\n")
		#print("recip_alleles2: "+str(recip_alleles2)+"\n")
		#print("donor_alleles1: "+str(donor_alleles1)+"\n")
		#print("donor_alleles2: "+str(donor_alleles2)+"\n")
		recip_allele2 = recip_alleles2[i-1]
		donor_allele1 = donor_alleles1[i-1]
		donor_allele2 = donor_alleles2[i-1]	

		for index,row in happair_selected.iterrows():
			mm_loc_count = 0
			mm_loc_sfvt_count = 0
			for pos in range(full_start, full_end):
				mm_aa_column_name = "MM_"+loc+"_"+str(pos)
				#happair_selected[mm_aa_column_name] = ''
				#print(donor_allele1[index],donor_allele2[index],recip_allele1[index],recip_allele2[index],pos)
				#allele1 = donor_allele1[index]
				#allele2 = donor_allele2[index]
				#allele3 = recip_allele1[index]
				#allele4 = recip_allele2[index]
				#print(type(allele1))
				#position = int(pos)
				#print(allele1)
				#print(type(position))
				#print(position)

				if (recip_allele1[index] == 'DRBX*NNNN' and recip_allele2[index] == 'DRBX*NNNN'):
					mm_pos_count == 2
				elif (donor_allele1[index] == 'DRBX*NNNN' and donor_allele2[index] == 'DRBX*NNNN' ):
					mm_pos_count == 2
				else:
					if (donor_allele1[index]  == 'DRBX*NNNN'):
						donor_allele1[index]  = donor_allele2[index] 
					if (donor_allele2[index]  == 'DRBX*NNNN'):
						donor_allele2[index]  = donor_allele1[index] 
					if (recip_allele1[index] == 'DRBX*NNNN'):
						recip_allele1[index] = recip_allele2[index] 
					if (recip_allele2[index]  == 'DRBX*NNNN'):
						recip_allele2[index]  = recip_allele1[index]
					mm_pos_count = genie.countAAMismatchesAllele(donor_allele1[index],donor_allele2[index],recip_allele1[index],recip_allele2[index],pos)
				mm_loc_count += mm_pos_count
				happair_selected.at[index,mm_aa_column_name] = mm_pos_count
				# happair_selected.at[index,mm_aa_column_name] = count_AA_Mismatches_Allele(HLA_seq,donor_allele1[index],donor_allele2[index],recip_allele1[index],recip_allele2[index],pos)

				# add columns for per-locus count
			mm_loc_count_name = "MM_" + loc + "_COUNT"
			happair_selected.at[index,mm_loc_count_name] = mm_loc_count
		
		'''
		mm_loc_count = 0
		mm_loc_sfvt_count = 0
		for pos in range(ard_start,ard_end):
			mm_aa_column_name = "MM_" + loc + "_" + str(pos)

			# initialize new columns
			happair_selected[mm_aa_column_name] = ''
			#happair_selected[mm_aa_column_name] = happair_selected["HAPPAIR_DONOR"]

			# get AA position and count mismatches for each row
			for index,row in happair_selected.iterrows():
				# if (loc=="DRB345"):
				# 	mm_pos_count=1
				# 	continue
				mm_pos_count = aa_mm.count_AA_Mismatches_Allele(donor_allele1[index],donor_allele2[index],recip_allele1[index],recip_allele2[index],pos)
				if loc == 'A':
					print('Locus:{}\t\tPosition:{}\t\tmm_pos_count:{}'.format(loc, str(pos), str(mm_pos_count)))
				mm_loc_count += mm_pos_count
				happair_selected.at[index,mm_aa_column_name] = mm_pos_count
				# happair_selected.at[index,mm_aa_column_name] = count_AA_Mismatches_Allele(HLA_seq,donor_allele1[index],donor_allele2[index],recip_allele1[index],recip_allele2[index],pos)

		# add columns for per-locus count
		mm_loc_count_name = "MM_" + loc + "_COUNT"
		happair_selected.at[index,mm_loc_count_name] = mm_loc_count
		'''

		# add SFVT columns
		# loop through all the SFVT_IDs
		for id in range(1,1000):
			seqf_name= "SFVT_" + loc + "_"  + str(id)
			if (seqf_name not in SFVT_list):
				continue
			# initialize new column for SFVT
			happair_selected[seqf_name] = happair_selected["HAPPAIR_DONOR"]			
			
			position_list = SFVT_positions[seqf_name]
			positions = [int(x) for x in position_list.split(",")]

			for index,row in happair_selected.iterrows():
				mm_pos_count = genie.countAAMismatchesAllele(donor_allele1[index],donor_allele2[index],recip_allele1[index],recip_allele2[index],positions)
				mm_loc_sfvt_count = mm_loc_count + mm_pos_count
				happair_selected.at[index,seqf_name] = mm_pos_count

			# add columns for per-locus count
			mm_loc_count_name = "MM_" + loc + "SFVT_COUNT"
			happair_selected.at[index,mm_loc_count_name] = mm_loc_sfvt_count

			# print ("Sequence Feature Name: " + seqf_name + " Positions: " + position_list + " MM Count: " + str(SFVT_MM_count))				

	# make sure data types are the same for merge
	tx_ki['PX_ID'] = tx_ki['PX_ID'].astype('int64')
	happair_selected['PX_ID'] = happair_selected['PX_ID'].astype('int64')

	'''
	# show complete HLA table of positions and mismatches
	print ("HapPair_Selected")
	print (happair_selected)

	print ('TX_KI')
	print (tx_ki)
	'''
	# merge in outcomes data from TX_KI file
	# tx_ki_design_matrix = pd.merge(tx_ki,happair_selected,how='inner',on='PX_ID')
	tx_ki_design_matrix = tx_ki.merge(happair_selected,how='right',left_on="PX_ID",right_on="PX_ID")

	'''
	print ("Design Matrix")
	with pd.option_context("display.max_rows",None,"display.max_columns",None):
		print (tx_ki_design_matrix)
	'''

	# write SRTR design matrix
	design_matrix_filename = pathloc + "/SRTR_AA_MM_9loc_matrix_genie_" + str(rep) + ".txt"
	tx_ki_design_matrix.to_csv(design_matrix_filename,index=False,sep="\t")

exit()