#!/usr/bin/env python
import pandas
import numpy
import re
import gzip
from aa_matching_msf import *

aa_mm = AAMatch(dbversion=3420)

# load imputation output data into a pair dictionary
# impute.srtr.*.csv.gz

pops_rollup = ['AFA','ASN','CAU','HIS','NAM','MLT']


		########Overall frequnecy already calculated for haplo pairs in 9loc impute data 
		# compute cumulative genotype frequency totals per subject
happair_id_total = {} # cumulative genotype frequency total
# compute cumulative genotype frequency totals per subject
happair_probs = {} # HLA probability distribution
happair_hla = {} # HLA haplotype pair distribution
PXID_list = {}

for pop in pops_rollup:
	impute_out_filename = "impute.srtr." + pop + ".csv.gz"
	impute_out_file = gzip.open(impute_out_filename, "rt")
	for line in impute_out_file:
		(subject_id,rank,hap1,hap2,freq) = line.split(',')
		if (subject_id == "subject	_id"): # skip header row
			continue

		happair_freq = 0
		happair_freq = float(freq)

		if subject_id not in happair_id_total:
			happair_id_total[subject_id] = 0
		happair_id_total[subject_id] += float(happair_freq)
	impute_out_file.close()

	impute_out_file = gzip.open(impute_out_filename, "rt")
	for line in impute_out_file:
		(subject_id,rank,hap1,hap2,freq) = line.split(',')
		if (subject_id == "subject_id"): # skip header row
			continue
		happair = hap1 + "+" + hap2

		happair_freq = 0
		happair_total = 0
		happair_freq = float(freq)

		if subject_id not in happair_probs:
			happair_probs[subject_id] = list()

		# renormalize to make sure freqs add up to zero
		happair_probs[subject_id].append(happair_freq / happair_id_total[subject_id])
		if subject_id not in happair_hla:
			happair_hla[subject_id] = list()
		happair_hla[subject_id].append(happair)
		PX_ID = subject_id[1:]
		# print (PX_ID)
		PXID_list[PX_ID] = 1


# load antigen mismatch equivalency tables from OPTN Policies
# a typing with a broad antigen is an antigen MM to an allele
# that is a split antigen

# Initialize an empty dictionary
allele_antigen_list = {}


# load broad/split antigen mappings
OPTN_antigen_broad_split_filename = "./OPTN_antigen_broad_split.txt"
OPTN_antigen_broad_split_file = open(OPTN_antigen_broad_split_filename, 'r')
OPTN_antigen_broad_split = {} # broad/split antigen mappings - split antigens are stored as a list
for line in OPTN_antigen_broad_split_file:
	line = line.rstrip()
	if line.startswith("OPTN_broad_antigen"): # skip header
		continue
	print (line)
	(broad_antigen,split_antigens) = line.split("\t")
	split_antigen_list = split_antigens.split("/")
	OPTN_antigen_broad_split[broad_antigen] = split_antigen_list
OPTN_antigen_broad_split_file.close()


# Path to the data file
OPTN_antigens_to_alleles_filename = "./OPTN_antigens_to_alleles_CPRA.txt"

# Open the file for reading
with open(OPTN_antigens_to_alleles_filename, 'r') as OPTN_antigens_to_alleles_file:
    # Iterate over each line in the file
    for row in OPTN_antigens_to_alleles_file:
        # Remove trailing whitespace
        row = row.rstrip()

        # Split the row into antigen and alleles
        antigen, alleles_gl = row.split("\t")

        # Exclude parental antigens from antigen list
        # Eliminate redundancy in alleles mapping from Bw4/Bw6 and ARD
        if (antigen in ["Bw4","Bw6","DRB1*14:01","DRB1*14:54"]):
            continue

        # Split the alleles string and iterate over each allele
        for allele in alleles_gl.split("/"):

            # alleles may be assigned to broad antigens in first pass
            if (allele in allele_antigen_list):
                cur_antigen = allele_antigen_list[allele]
                # allow for antigen assignment to be overwritten only allele has been assigned a broad
                if (cur_antigen in OPTN_antigen_broad_split):
                    allele_antigen_list[allele] = antigen
            else: # allele hasn't been assigned to antigen
                # Assign the antigen as the value for each allele key
                allele_antigen_list[allele] = antigen

# Print the resulting dictionary
# print(allele_antigen_list)

OPTN_antigens_to_alleles_file.close()


# count antigen mismatches, adjusting for donor homozygosity
def antigen_mm(locus,don_typ1,don_typ2,rec_typ1,rec_typ2):

	if (don_typ1 == "MISSING"):
		return "NA"
	if (rec_typ1 == "MISSING"):
		return "NA"


	donor_homoz = 0
	if (don_typ1 == don_typ2):
		donor_homoz = 1

	# allele-specific rollup to equivalent antigen
	if (locus in ["C","DQ"]):
		if (don_typ1 in allele_antigen_list): 
			don_typ1 = allele_antigen_list[don_typ1]
		if (don_typ2 in allele_antigen_list): 
			don_typ2 = allele_antigen_list[don_typ2]
		if (rec_typ1 in allele_antigen_list): 
			rec_typ1 = allele_antigen_list[rec_typ1]
		if (rec_typ2 in allele_antigen_list): 
			rec_typ2 = allele_antigen_list[rec_typ2]		

	if (locus == "DQA1"):
		(don_typ1_family,don_typ1_allele) = don_typ1.split(":")
		don_typ1 = don_typ1_family
		(don_typ2_family,don_typ2_allele) = don_typ2.split(":")
		don_typ2 = don_typ2_family
		(rec_typ1_family,rec_typ1_allele) = rec_typ1.split(":")
		rec_typ1 = rec_typ1_family
		(rec_typ2_family,rec_typ2_allele) = rec_typ2.split(":")
		rec_typ2 = rec_typ2_family

	if (locus == "DPA1"):
		(don_typ1_family,don_typ1_allele) = don_typ1.split(":")
		don_typ1 = don_typ1_family
		(don_typ2_family,don_typ2_allele) = don_typ2.split(":")
		don_typ2 = don_typ2_family
		(rec_typ1_family,rec_typ1_allele) = rec_typ1.split(":")
		rec_typ1 = rec_typ1_family
		(rec_typ2_family,rec_typ2_allele) = rec_typ2.split(":")
		rec_typ2 = rec_typ2_family

	mm_count = 0
	if ((don_typ1 != rec_typ1) & (don_typ1 != rec_typ2)):
		mm_count+=1

	if ((don_typ2 != rec_typ1) & (don_typ2 != rec_typ2)):
		mm_count+=1

	if ((mm_count == 2) & (donor_homoz == 1)):
		mm_count = 1

	return mm_count



# output donor imputation input file


impute_err_filename = "impute_failed_ids.csv"
impute_err_file = open(impute_err_filename, "w")
impute_err_file.write("Type,ID")


multiple_imputation_replicates = 10
for rep in range(1,multiple_imputation_replicates+1):

	antigen_mm_filename = "srtr_antigen_mm_" + str(rep) + ".csv"
	antigen_mm_file = open(antigen_mm_filename, "w")

	antigen_mm_file.write ("PX_ID,DON_C1,DON_C2,REC_C1,REC_C2,REC_C_MM_EQUIV_CUR,DON_DQ1,DON_DQ2,REC_DQ1,REC_DQ2,REC_DQ_MM_EQUIV_CUR,DON_DQA1,DON_DQA2,REC_DQA1,REC_DQA2,REC_DQA1_MM_EQUIV_CUR,DON_DPA1,DON_DPA2,REC_DPA1,REC_DPA2,REC_DPA1_MM_EQUIV_CUR\n")


	for PX_ID in PXID_list:

		SUBJECT_ID_RECIP = "R" + PX_ID
		SUBJECT_ID_DONOR = "D" + PX_ID
		# check for missing imputation output
		if (SUBJECT_ID_RECIP not in happair_hla):
			# print ("Missing Recip ID: " + PERS_ID)
			if (rep == 1):
				impute_err_file.write(SUBJECT_ID_RECIP)
			continue

		if (SUBJECT_ID_DONOR not in happair_hla):
			# print ("Missing Donor ID: " + DONOR_ID)
			if (rep == 1):
				impute_err_file.write(SUBJECT_ID_DONOR)
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


		# print (ID)
		# print ("C Antigens - DONOR: " + DON_C1 + " " + DON_C2 + " RECIP: " + REC_C1 + " " + REC_C2)
		REC_C_MM_EQUIV_CUR = antigen_mm("C",c1_donor,c2_donor,c1_recip,c2_recip)
		# print ("C Antigen MM: " + str(REC_C_MM_EQUIV_CUR))
		# print ("DQ Antigens - DONOR: " + DON_DQ1 + " " + DON_DQ2 + " RECIP: " + REC_DQ1 + " " + REC_DQ2)
		REC_DQ_MM_EQUIV_CUR = antigen_mm("DQ",dqb1_1_donor,dqb1_2_donor,dqb1_1_recip,dqb1_2_recip)
		# print ("DQ Antigen MM: " + str(REC_DQ_MM_EQUIV_CUR))
		REC_DQA1_MM_EQUIV_CUR = antigen_mm("DQA1",dqa1_1_donor,dqa1_2_donor,dqa1_1_recip,dqa1_2_recip)

		REC_DPA1_MM_EQUIV_CUR = antigen_mm("DPA1",dpa1_1_donor,dpa1_2_donor,dpa1_1_recip,dpa1_2_recip)


		antigen_mm_file.write(','.join([PX_ID,c1_donor,c2_donor,c1_recip,c2_recip,str(REC_C_MM_EQUIV_CUR),dqb1_1_donor,dqb1_2_donor,dqb1_1_recip,dqb1_2_recip,str(REC_DQ_MM_EQUIV_CUR),dqa1_1_donor,dqa1_2_donor,dqa1_1_recip,dqa1_2_recip,str(REC_DQA1_MM_EQUIV_CUR),dpa1_1_donor,dpa1_2_donor,dpa1_1_recip,dpa1_2_recip,str(REC_DPA1_MM_EQUIV_CUR)]) + "\n")