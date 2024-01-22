#!/usr/bin/env python

import gzip
import sys
import pandas as pd
from collections import defaultdict
import os
import timeit
from os import path
from aa_matching_msf import *

loci  = ["A","B","C","DRB1","DQB1"]

ard_start_pos = {
    "A" : 1,
    "B" : 1,
    "C" : 1,
    "DRB1" : 1,
    "DQA1" : 1,
    "DQB1" : 1,
    "DPA1" : 1,
    "DPB1" : 1,
}
ard_end_pos = {
    "A" : 182,
    "B" : 182,
    "C" : 182,
    "DRB1" : 94,
    "DQA1" : 87,
    "DQB1" : 94,
    "DPA1" : 84,
    "DPB1" : 92,
}

# Utility function to create dictionary
def multi_dict(K, type):
    if K == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: multi_dict(K-1, type))
 

# subroutine that returns dict of TRS values for a set of subjects
def TRS_locus(GF):
    TRS_subject = {}
    TRS_mean = 0
    for subject_id in GF:
        # print (subject_id)
        TRS = 0  # typing resolution score
        for GENO in GF[subject_id]:
            prob = GF[subject_id][GENO]
            # print (GENO)
            # print (prob)
            TRS = TRS + (prob * prob)
        # print ("TRS")
        TRS = round(TRS,6)
        # print (TRS)
        TRS_subject[subject_id] = TRS

    for subject_id in TRS_subject:
        TRS_mean += TRS_subject[subject_id] / ndonor

    # returns dict of TRS values for each subject
    return (TRS_subject,TRS_mean)


# filter out individuals not included in outcomes study
subject_ID_ethnicity_study = {}  # race/ethnicity of subjects included study - includes donor/recip status
grffail_filename = "SRTR_AA_MM_matrix_grffail_1.txt"  # from SRTR_AA_MM_matrix_grffail_replicates_2022-04-03.tar.gz
grffail_file = open(grffail_filename, "r")

for line in grffail_file:
    (PX_ID,CAN_RACE,REC_A_MM_EQUIV_CUR,REC_B_MM_EQUIV_CUR,REC_DR_MM_EQUIV_CUR,DON_RACE,*rest) = line.split('\t')
    if (PX_ID == "PX_ID"): # skip header row
        continue

    donor_ID = "D" + PX_ID
    recip_ID = "R" + PX_ID

    # rename multi-racial to MLT
    if (DON_RACE == "Multi-Racial"):
        DON_RACE = "MLT"
    if (CAN_RACE == "Multi-Racial"):
        CAN_RACE = "MLT"

    subject_ID_ethnicity_study[donor_ID] = DON_RACE
    subject_ID_ethnicity_study[recip_ID] = CAN_RACE
    #print (donor_ID + " " + DON_RACE)

grffail_file.close()

nsubjects = len(subject_ID_ethnicity_study)
print ("Number of subjects in grffail file: " + str(nsubjects))

# compute the number of subjects by category to get denominators for averages

nsubject_ethnicity = defaultdict(int)
nsubject_donor_ethnicity = defaultdict(int)
nsubject_recip_ethnicity = defaultdict(int)
nsubject_donor = 0
nsubject_recip = 0

for subject_id in subject_ID_ethnicity_study:

    # ethnicity
    subject_ethnicity = subject_ID_ethnicity_study[subject_id]
    nsubject_ethnicity[subject_ethnicity] += 1

    # donor vs recip
    if subject_id.startswith("D"):
        nsubject_donor += 1
        nsubject_donor_ethnicity[subject_ethnicity] += 1
    else:
        nsubject_recip += 1
        nsubject_recip_ethnicity[subject_ethnicity] += 1

ethnicity_list = list(nsubject_ethnicity.keys())

print (ethnicity_list)
print ("Number of donors: " + str(nsubject_donor))
print ("Number of recipients: " + str(nsubject_recip))
print ("Number of CAU donors: " + str(nsubject_donor_ethnicity["CAU"]))
print ("Number of CAU recips: " + str(nsubject_recip_ethnicity["CAU"]))


loci_selected = ["DRB1"]
# DRB1_AA_positions_selected = [13,28,30]
DRB1_AA_positions_selected = [9,10,11,12,13,26,28,30]
IDs_selected = ["R856948", "R867757", "R656410"]

# load SRTR HapLogic imputation output file and load probabilities for all subjects
impute_outfilename = "./srtr_impute.stage.gz"
impute_outfile = gzip.open(impute_outfilename, "rt")

# compute cumulative genotype frequency totals per subject
happair_id_total = {} # cumulative genotype frequency total
for line in impute_outfile:
    (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split('\t')
    if (cohort_id == "cohort_id"): # skip header row
        continue

    # if (subject_id not in IDs_selected):
    # 	continue

    if (subject_id not in subject_ID_ethnicity_study):
        continue

    happair_freq = 0
    if (hap1 == hap2):
        happair_freq = float(freq1) * float(freq2)

    else:
        happair_freq = 2 * float(freq1) * float(freq2)

    if subject_id not in happair_id_total:
        happair_id_total[subject_id] = 0
    happair_id_total[subject_id] += float(happair_freq)

    # print (subject_id + " " + hap1 + "+" + hap2 + " " + str(happair_id_total[subject_id]))

impute_outfile.close()



impute_outfile = gzip.open(impute_outfilename, "rt")

# compute probabilties for each haplotype pair
happair_probs = {} # HLA probability distribution
happair_hla = {} # HLA haplotype pair distribution
for line in impute_outfile:
    (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split('\t')
    if (cohort_id == "cohort_id"): # skip header row
        continue

    # if (subject_id not in IDs_selected):
    # 	continue

    if (subject_id not in subject_ID_ethnicity_study):
        continue

    happair = hap1 + "+" + hap2

    happair_freq = 0
    happair_total = 0
    if (hap1 == hap2):
        happair_freq = float(freq1) * float(freq2)

    else:
        happair_freq = 2 * float(freq1) * float(freq2)

    if subject_id not in happair_probs:
        happair_probs[subject_id] = list()
    happair_probs[subject_id].append(happair_freq / happair_id_total[subject_id])
    if subject_id not in happair_hla:
        happair_hla[subject_id] = list()
    happair_hla[subject_id].append(happair)
    # print (subject_id + " " + happair + " " + str(happair_freq))

# print (happair_hla)
# print (happair_probs)	

# create HLA Amino acid sequence database
aa_mm = AAMatch(dbversion=3420)

# compute probabilities for each AA-level genotype
TRS_dict = multi_dict(3, float) # TRS per subject per locus per position
subject_counter = 0
for subject_id in happair_hla:

    # starttime = timeit.default_timer()

    happair_index = 0
    AA_geno_probs = multi_dict(3, float)  # AA-level genotype list and probs per locus position
    for happair in happair_hla[subject_id]:
        # print (subject_id)
        # print ("Haplotype Pair: " + happair)
        (hap1,hap2) = happair.split('+')
        #print(hap1)
        (a1,c1,b1,drb1_1,dqb1_1) = hap1.split('~')
        (a2,c2,b2,drb1_2,dqb1_2) = hap2.split('~')

        if ((a1 == "A*11:52") | (a2 == "A*11:52")):
            a1 = "A*11:52Q"
            a2 = "A*11:52Q"

        if ((a1 == "A*23:19Q") | (a2 == "A*23:19Q")):
            a1 = "A*23:19N"
            a2 = "A*23:19N"

        if ((b1 == "B*13:08Q") | (b2 == "B*13:08Q")):
            b1 = "B*13:08"
            b2 = "B*13:08"

        if ((b1 == "B*15:22") | (b2 == "B*15:22")):
            b1 = "B*35:43"
            b2 = "B*35:43"

        if ((c1 == "C*15:20") | (c2 == "C*15:20")):
            c1 = "C*15:27"
            c2 = "C*15:27"

        if ((c1 == "C*03:12") | (c2 == "C*03:12")):
            c1 = "C*03:19"
            c2 = "C*03:19"

        if ((c1 == "C*03:23") | (c2 == "C*03:23")):
            c1 = "C*03:23N"
            c2 = "C*03:23N"

        if ((c1 == "C*05:02") | (c2 == "C*05:02")):
            c1 = "C*05:09"
            c2 = "C*05:09"

        if ((c1 == "C*13:01") | (c2 == "C*13:01")):
            c1 = "C*12:02"
            c2 = "C*12:02"

        if ((drb1_1 == "DRB1*07:02") | (drb1_2 == "DRB1*07:02")):
            drb1_1 = "DRB1*07:01"
            drb1_2 = "DRB1*07:01"


        happair_prob_list = happair_probs[subject_id]
        # print (happair_prob_list)
        happair_prob = happair_probs[subject_id][happair_index]
        happair_index += 1

        for locus in loci:
        # for locus in loci_selected:

            if (locus == "A"):
                allele1 = a1
                allele2 = a2
            if (locus == "B"):
                allele1 = b1
                allele2 = b2
            if (locus == "C"):
                allele1 = c1
                allele2 = c2
            if (locus == "DRB1"):
                allele1 = drb1_1
                allele2 = drb1_2
            if (locus == "DQB1"):
                allele1 = dqb1_1
                allele2 = dqb1_2

            if subject_id not in happair_probs:
                happair_probs[subject_id] = list()

            # get appropriate position range per locus
            # lots of incomplete sequences in IMGT/HLA that are filled in - use only ARD positions
            # for position in DRB1_AA_positions_selected:
            for position in range(ard_start_pos[locus],ard_end_pos[locus]):
                # print (locus)
                # print ("DRB1 allele-level genotypes: " + drb1_1 + "+" + drb1_2)

                AA1 = aa_mm.getAAposition(allele1,position)
                AA2 = aa_mm.getAAposition(allele2,position)
                (AA1,AA2) = sorted([AA1,AA2]) # alpha sort positions
                AA_geno = AA1 + "+" + AA2
                # print ("AA position: " + str(position))
                # print ("AA-level genotype: " + AA_geno)
                # print ("Hap Pair Prob: " + str(happair_prob))
                AA_geno_probs[locus][position][AA_geno] += happair_prob

    for locus in loci:
    # for locus in loci_selected:
        # print (locus)
        # for position in DRB1_AA_positions_selected:
        for position in range(ard_start_pos[locus],ard_end_pos[locus]):
            # print ("AA position: " + str(position))
            TRS = 0
            for AA_geno in AA_geno_probs[locus][position]:
                AA_geno_prob = AA_geno_probs[locus][position][AA_geno]
                AA_geno_prob = round(AA_geno_prob,10)
                # print ("Amino acid level genotype and probs: " + AA_geno + ": " + str(AA_geno_prob))
                TRS = TRS + (AA_geno_prob * AA_geno_prob)
                TRS_increment = AA_geno_prob * AA_geno_prob
                # print ("AA genotype contribution to TRS: " + str(TRS_increment))
            TRS_dict[subject_id][locus][position] = TRS

    # print("The time difference is :", timeit.default_timer() - starttime)
    subject_counter += 1
    if ((subject_counter % 10000) == 0):
        print ("Number of subjects with TRS computed: " + str(subject_counter))



# print out summary of all TRS values
TRS_average = multi_dict(4, float) # Average TRS per donor/recip per race per locus per position

# compute averages TRS per position for donors/recips
# compute averages TRS per position by race/ethnicity
# compute average TRS per position across the total multiethnic dataset
# HIS_count = 0
for subject_id in subject_ID_ethnicity_study:

    donor_or_recip = ""
    if (subject_id.startswith("D")):
        donor_or_recip = "DONOR"
    else:
        donor_or_recip = "RECIP"

    ethnicity = subject_ID_ethnicity_study[subject_id]

    # print (subject_id + " " + ethnicity + " " + donor_or_recip)


    # if (ethnicity == "HIS"):
    #     HIS_count += 1
    #     print ("Hispanic count: " +  str(HIS_count))

    for locus in loci:	
    # for locus in loci_selected:
        for position in range(ard_start_pos[locus],ard_end_pos[locus]):
        # for position in DRB1_AA_positions_selected:
            TRS = TRS_dict[subject_id][locus][position]
            # print (locus + " " + str(position) + " " + str(TRS))
            # print ("nsubjects: " + str(nsubjects))
            # print ("nsubjets of ethnicity " + ethnicity + " : " + str(nsubject_ethnicity[ethnicity]))
            # print ("ndonors: " + str(nsubject_donor))
            # print ("nrecips: " + str(nsubject_recip))
            # print ("ndonors of ethnicity " + ethnicity + " : " + str(nsubject_donor_ethnicity[ethnicity]))
            # print ("nrecips of ethnicity " + ethnicity + " : " + str(nsubject_recip_ethnicity[ethnicity]))
            TRS_average["SUBJECT"]["ALL"][locus][position] += (TRS / nsubjects)
            TRS_average["SUBJECT"][ethnicity][locus][position] += (TRS / nsubject_ethnicity[ethnicity])
            if (donor_or_recip == "DONOR"):
                TRS_average["DONOR"]["ALL"][locus][position] += (TRS / nsubject_donor)
                TRS_average["DONOR"][ethnicity][locus][position] += (TRS / nsubject_donor_ethnicity[ethnicity]) 
            else:
                TRS_average["RECIP"]["ALL"][locus][position] += (TRS / nsubject_recip)
                TRS_average["RECIP"][ethnicity][locus][position] += (TRS / nsubject_recip_ethnicity[ethnicity])
            # print ("Ethnicity: " + ethnicity)
            # print ("TRS Average for ALL Subjects Ethnicity: " + str(TRS_average["ALL_SUBJECTS"][ethnicity][locus][position]))
            # if (TRS_average["ALL_SUBJECTS"][ethnicity][locus][position] > 1):
            #     exit()

# print summary table
ethnicity_list.append("ALL")
TRS_average_out_filename = "SRTR_HLA_AA_TRS_Average.csv"
TRS_average_out_file = open(TRS_average_out_filename,"w")
TRS_average_out_file.write("Subject_Type,Ethnicity,Locus,AA_Position,TRS_Average\n")
for donor_or_recip in ["SUBJECT","DONOR","RECIP"]:
    for ethnicity in ethnicity_list:
        # for locus in loci_selected:
        for locus in loci:
            for position in range(ard_start_pos[locus],ard_end_pos[locus]):
            # for position in DRB1_AA_positions_selected:
                print (donor_or_recip + " " + ethnicity + " " + locus + " " + str(position))
            
                # UNK recip has 0 cases, so TRS should be NA
                TRS = ""
                if (not TRS_average[donor_or_recip][ethnicity][locus][position]):
                    TRS = "NA"
                else:
                    TRS = str(round(TRS_average[donor_or_recip][ethnicity][locus][position],8))
                
                print ("Average TRS: " + TRS)
                TRS_out_string = ",".join([donor_or_recip,ethnicity,locus,str(position),TRS])
                TRS_average_out_file.write(TRS_out_string + "\n")

# TODO - create separate table per locus
# rows in each tables are donor/recip-ethnicity combo
# columns are positions


exit()