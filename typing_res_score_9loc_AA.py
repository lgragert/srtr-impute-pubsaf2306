#!/usr/bin/env python

import gzip
from itertools import count
from pathlib import Path
import sys
import pandas as pd
from collections import defaultdict
import random
import pandas as pd
import re
import gzip
import numpy as np
import os
import glob
import pprint
import os
from time import time
from os import path, sep
from aa_matching_msf import *
aa_mm = AAMatch(dbversion=3420)


# Compute typing resolution score
#pathloc = str(path(__file__).parent) + '/'
#pathloc = sys.argv[1]
pathloc = str(Path('.'))
population_set = "ALL"

loci_value = 4
loci_range_value = 5

recip_alleles1 = []
recip_alleles2 = []
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

    # subroutine that returns dict of TRS values for a set of subjects
def DTRS_locus(DGF):
    DTRS_subject = {}
    DTRS_mean = 0
    for subject_id in DGF:
        # print (subject_id)
        DTRS = 0  # typing resolution score
        for DGENO in DGF[subject_id]:
            Dprob = DGF[subject_id][DGENO]
            # print (GENO)
            # print (prob)
            DTRS = DTRS + (Dprob * Dprob)
        # print ("TRS")
        DTRS = round(DTRS,6)
        # print (TRS)
        DTRS_subject[subject_id] = DTRS

    for subject_id in DTRS_subject:
        DTRS_mean += DTRS_subject[subject_id] / Dndonor

    # returns dict of TRS values for each subject
    return (DTRS_subject,DTRS_mean)


#TRS for AA
# Python strings start at zero, but amino acid sequences start at position 1
def getAAsubstring(hla_seq_dict,allele,start_position,end_position):
    return hla_seq_dict[allele].seq[start_position-1:end_position]


def getAAposition(hla_seq_dict,allele,position):
	return hla_seq_dict[allele].seq[position-1]


def TRS_AA(valuesList):
    TRS = 0  # typing resolution score
    for prob in valuesList:
        TRS = TRS + (prob * prob)
    TRS = round(TRS,6)
    return(TRS)

pops_NEMO = ['AFA','API','CAU','HIS','NAM','AAFA','AFB','CARB','SCAMB','AINDI','FILII','HAWI','JAPI','KORI','NCHI','SCSEAI','VIET','EURCAU','MENAFC','CARHIS','MSWHIS','SCAHIS','AISC','ALANAM','AMIND','CARIBI']
pops_US = ['AFA','API','CAU','HIS','NAM','DEC','MLT','OTH','UNK','AAFA','AFB','CARB','SCAMB','AINDI','FILII','HAWI','JAPI','KORI','NCHI','SCSEAI','VIET','EURCAU','MENAFC','CARHIS','MSWHIS','SCAHIS','AISC','ALANAM','AMIND','CARIBI']
pops_GLOBAL = ['US','DE','BR','IL','GB','PL','CN','TR','IN','ES','JP','CA','PT','NL','IT','FR','TH','AR','GR','SE','CY','CH','AU','SG','CL','CZ','SA','MX','AT','BE','HR','DK','RO','FI','IR','ZA','AM','NO','RU','SI','IE','SK','LT','NZ','RS','HU','BG','MK','UY','NG','GLOBAL']
pops_CPRA = ['AFA','CAU','HIS','NAM','ASN','HPI','MLT','HISETH']
pops_ALL = ['AFA','API','CAU','HIS','NAM']
patient = ['recip','donor']

trs = []
Dtrs=[]

#AA trs dict
trs_aa = pd.DataFrame()
Dtrs_aa = pd.DataFrame()
trs_aa_test = pd.DataFrame()


if (population_set == "ALL"):
    pops = pops_ALL
#print ("Pop,TRS_A,TRS_C,TRS_B,TRS_DRB1,TRS_DQB1")

for pop in pops:

    for DR in patient:

        print(DR)
        print (pop)
        impute_outfilename = pathloc + "/5locimpute_srtr_" + pop + "_" + DR + ".csv"
        print(impute_outfilename)
        impute_outfile = open(impute_outfilename, "rt")
        #with open(impute_outfilename, "rt") as myfile:
            #impute_outfile = [next(myfile) for x in range(100)]


        # Overall frequency already calculated for haplo pairs in 9loc impute data

        # compute cumulative genotype frequency totals per subject for later renormalization
        happair_id_total = {} # cumulative genotype frequency total
        Dhappair_id_total = {}

        aa_total = pd.DataFrame

        happair_selected_donor = {} # haplotype pair selected by weighted choice
        happair_selected_recip = {} # haplotype pair selected by weighted choice

        if (DR == "recip"):
            for line in impute_outfile:
                (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split(',')
                if (cohort_id == "cohort_id"): # skip header row
                    continue
                happair_freq = 0
                if (hap1 == hap2):
                    happair_freq = float(freq1) * float(freq2)

                else:
                    happair_freq = 2 * float(freq1) * float(freq2)

                if subject_id not in happair_id_total:
                    happair_id_total[subject_id] = 0
                happair_id_total[subject_id] += float(happair_freq)

            impute_outfile.close()
            
            ndonor = len(happair_id_total)

            # reopen impute file for 2nd pass
            impute_outfile = open(impute_outfilename, "rt")

            # compute probabilties for each haplotype pair
            happair_probs = {} # HLA probability distribution
            happair_hla = {} # HLA haplotype pair distribution
            GF_A = defaultdict(dict)
            GF_B = defaultdict(dict)
            GF_C = defaultdict(dict)
            GF_DRB1 = defaultdict(dict)
            GF_DQB1 = defaultdict(dict)
            AA_C = defaultdict(dict)
            AA_B = defaultdict(dict)
            AA_DRB1 = defaultdict(dict)
            AA_DQB1 = defaultdict(dict)
            AA_A = defaultdict(dict)
            AA_C = defaultdict(dict)
            AA_B = defaultdict(dict)
            AA_DRB1 = defaultdict(dict)
            AA_DQB1 = defaultdict(dict)

            for line in impute_outfile:
                (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split(',')
                if (cohort_id == "cohort_id"): # skip header row
                    continue
                (A_1,C_1,B_1,DRB1_1,DQB1_1) = hap1.split ('~')
                (A_2,C_2,B_2,DRB1_2,DQB1_2) = hap2.split ('~')        

                happair = hap1 + "+" + hap2
                happair_freq = 0
                if (hap1 == hap2):
                    happair_freq = float(freq1) * float(freq2)

                else:
                    happair_freq = 2 * float(freq1) * float(freq2)
                    
                # sort strings
                (A_1,A_2) = sorted ([A_1,A_2])
                (B_1,B_2) = sorted ([B_1,B_2])
                (C_1,C_2) = sorted ([C_1,C_2])
                (DRB1_1,DRB1_2) = sorted ([DRB1_1,DRB1_2])
                (DQB1_1,DQB1_2) = sorted ([DQB1_1,DQB1_2])

                GENO_A = '+'.join([A_1,A_2])
                GENO_C = '+'.join([C_1,C_2])
                GENO_B = '+'.join([B_1,B_2])
                GENO_DRB1 = '+'.join([DRB1_1,DRB1_2])
                GENO_DQB1 = '+'.join([DQB1_1,DQB1_2])

                # print (happair)
                # print (GENO_A)
                #print(happair_freq)
                #print(happair_id_total)
                #print(subject_id)
                
                prob = happair_freq / happair_id_total[subject_id]

                #print(prob)

                # print (prob)

                if GENO_A not in GF_A[subject_id]:
                    GF_A[subject_id][GENO_A] = prob
                else:
                    GF_A[subject_id][GENO_A] = GF_A[subject_id][GENO_A] + prob

                if GENO_B not in GF_B[subject_id]:
                    GF_B[subject_id][GENO_B] = prob
                else:
                    GF_B[subject_id][GENO_B] = GF_B[subject_id][GENO_B] + prob

                if GENO_C not in GF_C[subject_id]:
                    GF_C[subject_id][GENO_C] = prob
                else:
                    GF_C[subject_id][GENO_C] = GF_C[subject_id][GENO_C] + prob

                if GENO_DRB1 not in GF_DRB1[subject_id]:
                    GF_DRB1[subject_id][GENO_DRB1] = prob
                else:
                    GF_DRB1[subject_id][GENO_DRB1] = GF_DRB1[subject_id][GENO_DRB1] + prob

                if GENO_DQB1 not in GF_DQB1[subject_id]:
                    GF_DQB1[subject_id][GENO_DQB1] = prob
                else:
                    GF_DQB1[subject_id][GENO_DQB1] = GF_DQB1[subject_id][GENO_DQB1] + prob

                #print(A_trs)
            (TRS_A_subject,TRS_A_mean) = TRS_locus(GF_A)
            (TRS_B_subject,TRS_B_mean) = TRS_locus(GF_B)
            (TRS_C_subject,TRS_C_mean) = TRS_locus(GF_C)
            (TRS_DRB1_subject,TRS_DRB1_mean) = TRS_locus(GF_DRB1)
            (TRS_DQB1_subject,TRS_DQB1_mean) = TRS_locus(GF_DQB1)

            #print(TRS_A_subject)
            trs.append({
                'Population': pop,
                'TRS_A_mean': TRS_A_mean,
                'TRS_C_mean': TRS_C_mean,
                'TRS_B_mean': TRS_B_mean,
                'TRS_DRB1_mean': TRS_DRB1_mean,
                'TRS_DQB1_mean': TRS_DQB1_mean
            })
            trs_dataframe = pd.DataFrame(trs)
            design_matrix_filename = pathloc + "/trs_5loc_matrix_recip.csv"
            trs_dataframe.to_csv(design_matrix_filename,index=False,sep=",")
            impute_outfile.close()

            impute_outfile = open(impute_outfilename, "rt")

#for AA TRS Values
            for line in impute_outfile:
                (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split(',')
                if (cohort_id == "cohort_id"): # skip header row
                    continue


                #print(subject_id)

                if subject_id not in AA_A:
                    A_sub =GF_A.get(subject_id)
                    #print(A_sub)
                    #keysList = list(A_sub.keys())
                    #print(keysList)
                    valuesList_A = list(A_sub.values())
                    #print(valuesList_A)
                    for keys in A_sub:
                        #print(keys)
                        A_list = keys.split('+')
                        A_A = defaultdict(dict)
                        #print(A_list)
                        for A in A_list:
                            if (A == "A*23:19Q"):
                                A = "A*23:19N"
                            if (A == "A*11:52"):
                                A = "A*11:52Q"
                            for i in range(1,loci_range_value):
                                Aloc = GENO_A[0].split('*')[0]
                                #print(Aloc_1)
                                #print(loc)
                                A_ard_start = aa_mm.ard_start_pos[Aloc]
                                A_ard_end = aa_mm.ard_end_pos[Aloc]
                                aa_A1_full = aa_mm.getAAsubstring(A,A_ard_start,A_ard_end)
                                aa_A1= list(aa_A1_full)
                                for pos in range(A_ard_start,A_ard_end):
                                    SAAP_A = str(Aloc+"_"+str(pos))
                                    #aas_A1 = aa_mm.getAAposition(Aloc_1,pos)
                                    #print(type(aa_A))
                                    AA1 = aa_A1[pos-1]
                                    A_A[A][SAAP_A]=AA1
                                    #print(AA_A)
                                #end = time()
                                #total = end - start
                                #print("Time to assign AAs for one Subject ID:", total)
                    A_A = pd.DataFrame.from_dict(A_A, orient = 'index')
                    #print(A_A)
                    #for column in A_A:
                        #print(column)
                        #An = A_A.count(column)
                        #An = dict(An)

                        #print(n)
                        #end = time()
                        #total = end - start
                        #print("Time to find unique number of AAs for one Subject ID:", total)
                    #start = time()
                    for pos in range(A_ard_start,A_ard_end):
                        SAAP_A = str(Aloc+"_"+str(pos))
                        #for SAAP_A in An:
                            #print(SAAP_A)
                        Atrsmatch = A_A[SAAP_A].nunique()
                        if (Atrsmatch == 1):
                            prob = 1
                            #print("Prob:",prob)
                            AA_A[subject_id][SAAP_A]=prob
                        else:
                            prob = TRS_AA(valuesList_A)
                            #print("Prob:",prob)
                            AA_A[subject_id][SAAP_A]=prob
                                #print(AA_A)
                        #end = time()
                        #total = end - start
                    #end = time()
                    #time = end-start
                    #print("Time to find one unique values for on loci of a haplotype:", time)
                        #print("Time to assign TRS for each position for one Subject ID:", total)

                    C_sub =GF_C.get(subject_id)
                    valuesList_C = list(C_sub.values())
                    for keys in C_sub:
                        C_list = keys.split('+')
                        A_C = defaultdict(dict)
                        for C in C_list:
                            if (C == "C*15:20"):
                                C = "C*15:27"
                            if (C == "C*05:02"):
                                C = "C*05:09"
                            if (C == "C*03:12"):
                                C = "C*03:19"
                            if (C == "C*03:23"):
                                C = "C*03:23N"
                            if (C == "C*03:23"):
                                C = "C*03:23N"
                            if (C == "C*13:01"):
                                C = "C*12:02"
                            for i in range(1,loci_range_value):
                                Cloc = GENO_C[0].split('*')[0]
                                C_ard_start = aa_mm.ard_start_pos[Cloc]
                                C_ard_end = aa_mm.ard_end_pos[Cloc]
                                aa_C_full = aa_mm.getAAsubstring(C,C_ard_start,C_ard_end)
                                aa_C= list(aa_C_full)
                                for pos in range(C_ard_start,C_ard_end):
                                    SAAP_C = str(Cloc+"_"+str(pos))
                                    AC = aa_C[pos-1]
                                    A_C[C][SAAP_C]=AC
                    A_C = pd.DataFrame.from_dict(A_C, orient = 'index')
                    #for column in A_C:
                        #Cn = A_C.nunique(axis=0)
                        #Cn = dict(Cn)
                    for pos in range(C_ard_start,C_ard_end):
                        SAAP_C = str(Cloc+"_"+str(pos))
                        #for SAAP_C in Cn:
                        Ctrsmatch = A_C[SAAP_C].nunique()
                        if (Ctrsmatch == 1):
                            prob = 1
                            AA_C[subject_id][SAAP_C]=prob
                        else:
                            prob = TRS_AA(valuesList_C)
                            AA_C[subject_id][SAAP_C]=prob


                    B_sub =GF_B.get(subject_id)
                    valuesList_B = list(B_sub.values())
                    for keys in B_sub:
                        B_list = keys.split('+')
                        A_B = defaultdict(dict)
                        for B in B_list:
                            if (B == "B*15:22"):
                                B = "B*35:43"
                            if (B == "B*13:08Q"):
                                B = "B*13:08"                                
                            for i in range(1,loci_range_value):
                                Bloc = GENO_B[0].split('*')[0]
                                B_ard_start = aa_mm.ard_start_pos[Bloc]
                                B_ard_end = aa_mm.ard_end_pos[Bloc]
                                aa_B_full = aa_mm.getAAsubstring(B,B_ard_start,B_ard_end)
                                aa_B= list(aa_B_full)
                                for pos in range(B_ard_start,B_ard_end):
                                    SAAP_B = str(Bloc+"_"+str(pos))
                                    AB = aa_B[pos-1]
                                    A_B[B][SAAP_B]=AB
                    A_B = pd.DataFrame.from_dict(A_B, orient = 'index')
                    #for column in A_B:
                        #Bn = A_B.nunique(axis=0)
                        #Bn = dict(Bn)
                    for pos in range(B_ard_start,B_ard_end):
                        SAAP_B = str(Bloc+"_"+str(pos))
                        #for SAAP_B in Bn:
                            #Btrsmatch = Bn.get(SAAP_B)
                        Btrsmatch = A_B[SAAP_B].nunique()
                        if (Btrsmatch == 1):
                            prob = 1
                            AA_B[subject_id][SAAP_B]=prob
                        else:
                            prob = TRS_AA(valuesList_B)
                            AA_B[subject_id][SAAP_B]=prob

                    DRB1_sub =GF_DRB1.get(subject_id)
                    valuesList_DRB1 = list(DRB1_sub.values())
                    for keys in DRB1_sub:
                        DRB1_list = keys.split('+')
                        A_DRB1 = defaultdict(dict)
                        for DRB1 in DRB1_list:
                            # key error due to not being in IMGT/HLA
                            if (DRB1 == "DRB1*07:02"):
                                DRB1 = "DRB1*07:01"
                            for i in range(1,loci_range_value):
                                DRB1loc = GENO_DRB1[0:4].split('*')[0]
                                DRB1_ard_start = aa_mm.ard_start_pos[DRB1loc]
                                DRB1_ard_end = aa_mm.ard_end_pos[DRB1loc]
                                aa_DRB1_full = aa_mm.getAAsubstring(DRB1,DRB1_ard_start,DRB1_ard_end)
                                aa_DRB1= list(aa_DRB1_full)
                                for pos in range(DRB1_ard_start,DRB1_ard_end):
                                    SAAP_DRB1 = str(DRB1loc+"_"+str(pos))
                                    ADRB1 = aa_DRB1[pos-1]
                                    A_DRB1[DRB1][SAAP_DRB1]=ADRB1
                    A_DRB1 = pd.DataFrame.from_dict(A_DRB1, orient = 'index')
                    #for column in A_DRB1:
                        #DRB1n = A_DRB1.nunique(axis=0)
                        #DRB1n = dict(DRB1n)
                    for pos in range(DRB1_ard_start,DRB1_ard_end):
                        SAAP_DRB1 = str(DRB1loc+"_"+str(pos))
                            #for SAAP_DRB1 in DRB1n:
                                #DRB1trsmatch = DRB1n.get(SAAP_DRB1)
                        DRB1trsmatch = A_DRB1[SAAP_DRB1].nunique()
                        if (DRB1trsmatch == 1):
                            prob = 1
                            AA_DRB1[subject_id][SAAP_DRB1]=prob
                        else:
                            prob = TRS_AA(valuesList_DRB1)
                            AA_DRB1[subject_id][SAAP_DRB1]=prob

                    DQB1_sub =GF_DQB1.get(subject_id)
                    valuesList_DQB1 = list(DQB1_sub.values())
                    for keys in DQB1_sub:
                        DQB1_list = keys.split('+')
                        A_DQB1 = defaultdict(dict)
                        for DQB1 in DQB1_list:
                            for i in range(1,loci_range_value):
                                DQB1loc = GENO_DQB1[0:4].split('*')[0]
                                DQB1_ard_start = aa_mm.ard_start_pos[DQB1loc]
                                DQB1_ard_end = aa_mm.ard_end_pos[DQB1loc]
                                aa_DQB1_full = aa_mm.getAAsubstring(DQB1,DQB1_ard_start,DQB1_ard_end)
                                aa_DQB1= list(aa_DQB1_full)
                                for pos in range(DQB1_ard_start,DQB1_ard_end):
                                    SAAP_DQB1 = str(DQB1loc+"_"+str(pos))
                                    ADQB1 = aa_DQB1[pos-1]
                                    A_DQB1[DQB1][SAAP_DQB1]=ADQB1
                    A_DQB1 = pd.DataFrame.from_dict(A_DQB1, orient = 'index')
                    #for column in A_DQB1:
                        #DQB1n = A_DQB1.nunique(axis=0)
                        #DQB1n = dict(DQB1n)
                    for pos in range(DQB1_ard_start,DQB1_ard_end):
                        SAAP_DQB1 = str(DQB1loc+"_"+str(pos))
                        DQB1trsmatch = A_DQB1[SAAP_DQB1].nunique()
                            #for SAAP_DQB1 in DQB1n:
                                #DQB1trsmatch = DQB1n.get(SAAP_DQB1)
                        if (DQB1trsmatch == 1):
                            prob = 1
                            AA_DQB1[subject_id][SAAP_DQB1]=prob
                        else:
                            prob = TRS_AA(valuesList_DQB1)
                            AA_DQB1[subject_id][SAAP_DQB1]=prob
                else:
                    continue
            AA_A = pd.DataFrame.from_dict(AA_A, orient = 'index')
            #print(AA_A)
            #AA_dataframe = pd.DataFrame(AA_A)
            #AA_design_matrix_filename = pathloc +"/AAgeno_trs_5loc_aa_matrix_recip.csv"
            #AA_dataframe.to_csv(AA_design_matrix_filename,index= True,sep=",")
            AA_A=AA_A.mean(axis='index')
            trs_aa = trs_aa.append(pd.DataFrame([AA_A], index=[pop] ))
            #trs_aa =trs_aa.concat([AA_A],index=[pop])

            AA_C = pd.DataFrame.from_dict(AA_C, orient = 'index')
            AA_C=AA_C.mean(axis='index')
            trs_aa = trs_aa.append(pd.DataFrame([AA_C], index=[pop] ))
            #trs_aa =trs_aa.concat([AA_C],index=[pop])

            AA_B = pd.DataFrame.from_dict(AA_B, orient = 'index')
            AA_B=AA_B.mean(axis='index')
            trs_aa = trs_aa.append(pd.DataFrame([AA_B], index=[pop] ))
            #trs_aa =trs_aa.concat([AA_B],index=[pop])

            AA_DRB1 = pd.DataFrame.from_dict(AA_DRB1, orient = 'index')
            AA_DRB1=AA_DRB1.mean(axis='index')
            trs_aa = trs_aa.append(pd.DataFrame([AA_DRB1], index=[pop] ))
            #trs_aa =pd.concat([AA_DRB1],index=[pop])


            AA_DQB1 = pd.DataFrame.from_dict(AA_DQB1, orient = 'index')
            AA_DQB1=AA_DQB1.mean(axis='index')
            trs_aa = trs_aa.append(pd.DataFrame([AA_DQB1], index=[pop] ))
            #print(trs_aa)
            #trs_aa =trs_aa.concat([AA_DRB1],index=[pop])


            trs_dataframe = pd.DataFrame(trs_aa)
            design_matrix_filename = pathloc +"/trs_5loc_aa_matrix_recip.csv"
            trs_dataframe.to_csv(design_matrix_filename,index= True,sep=",")




        if (DR == 'donor'):
            for line in impute_outfile:
                (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split(',')
                if (cohort_id == "cohort_id"): # skip header row
                    continue
            
                Dhappair_freq = 0
                if (hap1 == hap2):
                    Dhappair_freq = float(freq1) * float(freq2)

                else:
                    Dhappair_freq = 2 * float(freq1) * float(freq2)

                if subject_id not in Dhappair_id_total:
                    Dhappair_id_total[subject_id] = 0
                Dhappair_id_total[subject_id] += float(Dhappair_freq)

            impute_outfile.close()
        
            Dndonor = len(Dhappair_id_total)

            # reopen impute file for 2nd pass
            impute_outfile = open(impute_outfilename, "rt")

            # compute probabilties for each haplotype pair
            Dhappair_probs = {} # HLA probability distribution
            Dhappair_hla = {} # HLA haplotype pair distribution
            DGF_A = defaultdict(dict)
            DGF_B = defaultdict(dict)
            DGF_C = defaultdict(dict)
            DGF_DRB1 = defaultdict(dict)
            DGF_DQB1 = defaultdict(dict)
            DAA_A = defaultdict(dict)
            DAA_C = defaultdict(dict)
            DAA_B = defaultdict(dict)
            DAA_DRB1 = defaultdict(dict)
            DAA_DQB1 = defaultdict(dict)
            for line in impute_outfile:
                (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split(',')
                if (cohort_id == "cohort_id"): # skip header row
                    continue

                (A_1,C_1,B_1,DRB1_1,DQB1_1) = hap1.split ('~')
                (A_2,C_2,B_2,DRB1_2,DQB1_2) = hap2.split ('~')        

                Dhappair = hap1 + "+" + hap2
                Dhappair_freq = 0
                if (hap1 == hap2):
                    Dhappair_freq = float(freq1) * float(freq2)

                else:
                    Dhappair_freq = 2 * float(freq1) * float(freq2)


                # sort strings
                (A_1,A_2) = sorted ([A_1,A_2])
                (B_1,B_2) = sorted ([B_1,B_2])
                (C_1,C_2) = sorted ([C_1,C_2])
                (DRB1_1,DRB1_2) = sorted ([DRB1_1,DRB1_2])
                (DQB1_1,DQB1_2) = sorted ([DQB1_1,DQB1_2])

                DGENO_A = '+'.join([A_1,A_2])
                DGENO_C = '+'.join([C_1,C_2])
                DGENO_B = '+'.join([B_1,B_2])
                DGENO_DRB1 = '+'.join([DRB1_1,DRB1_2])
                DGENO_DQB1 = '+'.join([DQB1_1,DQB1_2])

                # print (happair)
                # print (GENO_A)
                
                Dprob = Dhappair_freq / Dhappair_id_total[subject_id]

                # print (prob)

                if DGENO_A not in DGF_A[subject_id]:
                    DGF_A[subject_id][DGENO_A] = Dprob
                else:
                    DGF_A[subject_id][DGENO_A] = DGF_A[subject_id][DGENO_A] + Dprob

                if DGENO_B not in DGF_B[subject_id]:
                    DGF_B[subject_id][DGENO_B] = Dprob
                else:
                    DGF_B[subject_id][DGENO_B] = DGF_B[subject_id][DGENO_B] + Dprob

                if DGENO_C not in DGF_C[subject_id]:
                    DGF_C[subject_id][DGENO_C] = Dprob
                else:
                    DGF_C[subject_id][DGENO_C] = DGF_C[subject_id][DGENO_C] + Dprob

                if DGENO_DRB1 not in DGF_DRB1[subject_id]:
                    DGF_DRB1[subject_id][DGENO_DRB1] = Dprob
                else:
                    DGF_DRB1[subject_id][DGENO_DRB1] = DGF_DRB1[subject_id][DGENO_DRB1] + Dprob

                if DGENO_DQB1 not in DGF_DQB1[subject_id]:
                    DGF_DQB1[subject_id][DGENO_DQB1] = Dprob
                else:
                    DGF_DQB1[subject_id][DGENO_DQB1] = DGF_DQB1[subject_id][DGENO_DQB1] + Dprob


            (DTRS_A_subject,DTRS_A_mean) = DTRS_locus(DGF_A)
            (DTRS_B_subject,DTRS_B_mean) = DTRS_locus(DGF_B)
            (DTRS_C_subject,DTRS_C_mean) = DTRS_locus(DGF_C)
            (DTRS_DRB1_subject,DTRS_DRB1_mean) = DTRS_locus(DGF_DRB1)
            (DTRS_DQB1_subject,DTRS_DQB1_mean) = DTRS_locus(DGF_DQB1)

            #print("DONOR:", dict(list(Dhappair_id_total.items())[0:10]), dict(list(DGF_A.items())[0:10]), dict(list(Dhappair_hla.items())[0:10]))


            # print (TRS_A_subject)
            # print ("TRS_A_mean")
            # print (TRS_A_mean)
            Dtrs.append({
                'Population': pop,
                'TRS_A_mean': DTRS_A_mean,
                'TRS_C_mean': DTRS_C_mean,
                'TRS_B_mean': DTRS_B_mean,
                'TRS_DRB1_mean': DTRS_DRB1_mean,
                'TRS_DQB1_mean': DTRS_DQB1_mean,
            })
            Dtrs_dataframe = pd.DataFrame(Dtrs)
            Ddesign_matrix_filename = pathloc + "/trs_5loc_matrix_donor.csv"
            Dtrs_dataframe.to_csv(Ddesign_matrix_filename,index=False,sep=",")

            impute_outfile.close()

            impute_outfile = open(impute_outfilename, "rt")

#for AA TRS Values
            for line in impute_outfile:
                (cohort_id,subject_id,pop1,hap1,freq1,pop2,hap2,freq2,freq_below_cutoff) = line.split(',')
                if (cohort_id == "cohort_id"): # skip header row
                    continue

                if subject_id not in DAA_A:
                    DA_sub =DGF_A.get(subject_id)
                    #print(A_sub)
                    #keysList = list(A_sub.keys())
                    #print(keysList)
                    valuesList_DA = list(DA_sub.values())
                    #print(valuesList_A)
                    for keys in DA_sub:
                        #print(keys)
                        DA_list = keys.split('+')
                        DA_A = defaultdict(dict)
                        #print(A_list)
                        for DA in DA_list:
                            if (DA == "A*23:19Q"):
                                DA = "A*23:19N"
                            if (DA == "A*11:52"):
                                DA = "A*11:52Q"
                            for i in range(1,loci_range_value):
                                DAloc = DGENO_A[0].split('*')[0]
                                #print(Aloc_1)
                                #print(loc)
                                DA_ard_start = aa_mm.ard_start_pos[DAloc]
                                DA_ard_end = aa_mm.ard_end_pos[DAloc]
                                Daa_A1_full = aa_mm.getAAsubstring(DA,DA_ard_start,DA_ard_end)
                                Daa_A1= list(Daa_A1_full)
                                for pos in range(DA_ard_start,DA_ard_end):
                                    SAAP_A = str(Aloc+"_"+str(pos))
                                    #aas_A1 = aa_mm.getAAposition(Aloc_1,pos)
                                    #print(type(aa_A))
                                    DAA1 = Daa_A1[pos-1]
                                    DA_A[DA][SAAP_A]=DAA1
                                        #print(AA_A)
                    DA_A = pd.DataFrame.from_dict(DA_A, orient = 'index')
                    #print(A_A)
                    #for column in A_A:
                        #print(column)
                        #DAn = DA_A.nunique(axis=0)
                        #DAn = dict(DAn)
                        #print(n)
                    for pos in range(DA_ard_start,DA_ard_end):
                        SAAP_A = str(DAloc+"_"+str(pos))
                        DAtrsmatch = DA_A[SAAP_A].nunique()

                            #for SAAP_A in DAn:
                                #print(SAAP_A)
                                #DAtrsmatch =DAn.get(SAAP_A)
                                #print(trsmatch)
                        if (DAtrsmatch == 1):
                            prob = 1
                            DAA_A[subject_id][SAAP_A]=prob
                        else:
                            prob = TRS_AA(valuesList_DA)
                            #print(prob)
                            DAA_A[subject_id][SAAP_A]=prob
                            #print(AA_A)

                    DC_sub =DGF_C.get(subject_id)
                    valuesList_DC = list(DC_sub.values())
                    for keys in DC_sub:
                        DC_list = keys.split('+')
                        DA_C = defaultdict(dict)
                        for DC in DC_list:
                            if (DC == "C*15:20"):
                                DC = "C*15:27"
                            if (DC == "C*05:02"):
                                DC = "C*05:09"
                            if (DC == "C*03:12"):
                                DC = "C*03:19"
                            if (DC == "C*03:23"):
                                DC = "C*03:23N"
                            if (DC == "C*03:23"):
                                DC = "C*03:23N"
                            if (DC == "C*13:01"):
                                DC = "C*12:02"
                            for i in range(1,loci_range_value):
                                DCloc = DGENO_C[0].split('*')[0]
                                DC_ard_start = aa_mm.ard_start_pos[DCloc]
                                DC_ard_end = aa_mm.ard_end_pos[DCloc]
                                aa_DC_full = aa_mm.getAAsubstring(DC,DC_ard_start,DC_ard_end)
                                aa_DC= list(aa_DC_full)
                                for pos in range(DC_ard_start,DC_ard_end):
                                    SAAP_C = str(DCloc+"_"+str(pos))
                                    DAC = aa_DC[pos-1]
                                    DA_C[DC][SAAP_C]=DAC
                    DA_C = pd.DataFrame.from_dict(DA_C, orient = 'index')
                    #for column in DA_C:
                        #DCn = DA_C.nunique(axis=0)
                        #DCn = dict(DCn)
                    for pos in range(C_ard_start,C_ard_end):
                        SAAP_C = str(DCloc+"_"+str(pos))
                        DCtrsmatch = DA_C[SAAP_C].nunique()
                            #for SAAP_C in DCn:
                                #DCtrsmatch =DCn.get(SAAP_C)
                        if (DCtrsmatch == 1):
                            prob = 1
                            DAA_C[subject_id][SAAP_C]=prob
                        else:
                            prob = TRS_AA(valuesList_DC)
                            DAA_C[subject_id][SAAP_C]=prob


                    DB_sub =DGF_B.get(subject_id)
                    valuesList_DB = list(DB_sub.values())
                    for keys in DB_sub:
                        DB_list = keys.split('+')
                        DA_B = defaultdict(dict)
                        for DB in DB_list:
                            if (DB == "B*15:22"):
                                DB = "B*35:43"
                            if (DB == "B*13:08Q"):
                                DB = "B*13:08"   
                        for i in range(1,loci_range_value):
                            DBloc = DGENO_B[0].split('*')[0]
                            DB_ard_start = aa_mm.ard_start_pos[DBloc]
                            DB_ard_end = aa_mm.ard_end_pos[DBloc]
                            aa_DB_full = aa_mm.getAAsubstring(DB,DB_ard_start,DB_ard_end)
                            aa_DB= list(aa_DB_full)
                            for pos in range(DB_ard_start,DB_ard_end):
                                SAAP_B = str(DBloc+"_"+str(pos))
                                DAB = aa_DB[pos-1]
                                DA_B[DB][SAAP_B]=DAB
                    DA_B = pd.DataFrame.from_dict(DA_B, orient = 'index')
                    #for column in A_B:
                        #DBn = DA_B.nunique(axis=0)
                        #DBn = dict(DBn)
                    for pos in range(DB_ard_start,DB_ard_end):
                        SAAP_B = str(DBloc+"_"+str(pos))
                        DBtrsmatch = DA_B[SAAP_B].nunique()
                            #for SAAP_B in DBn:
                                #DBtrsmatch =DBn.get(SAAP_B)
                        if (DBtrsmatch == 1):
                            prob = 1
                            DAA_B[subject_id][SAAP_B]=prob
                        else:
                            prob = TRS_AA(valuesList_DB)
                            DAA_B[subject_id][SAAP_B]=prob

                    DDRB1_sub =DGF_DRB1.get(subject_id)
                    valuesList_DDRB1 = list(DDRB1_sub.values())
                    for keys in DDRB1_sub:
                        DDRB1_list = keys.split('+')
                        DA_DRB1 = defaultdict(dict)
                        for DDRB1 in DDRB1_list:
                            # key error due to not being in IMGT/HLA
                            if (DDRB1 == "DRB1*07:02"):
                                DDRB1 = "DRB1*07:01"
                            for i in range(1,loci_range_value):
                                DDRB1loc = DGENO_DRB1[0:4].split('*')[0]
                                DDRB1_ard_start = aa_mm.ard_start_pos[DDRB1loc]
                                DDRB1_ard_end = aa_mm.ard_end_pos[DDRB1loc]
                                aa_DDRB1_full = aa_mm.getAAsubstring(DDRB1,DDRB1_ard_start,DDRB1_ard_end)
                                aa_DDRB1= list(aa_DDRB1_full)
                                for pos in range(DDRB1_ard_start,DDRB1_ard_end):
                                    SAAP_DRB1 = str(DDRB1loc+"_"+str(pos))
                                    DADRB1 = aa_DDRB1[pos-1]
                                    DA_DRB1[DDRB1][SAAP_DRB1]=DADRB1
                    DA_DRB1 = pd.DataFrame.from_dict(DA_DRB1, orient = 'index')
                    #for column in DA_DRB1:
                        #DDRB1n = DA_DRB1.nunique(axis=0)
                        #DDRB1n = dict(DDRB1n)
                    for pos in range(DDRB1_ard_start,DDRB1_ard_end):
                        SAAP_DRB1 = str(DDRB1loc+"_"+str(pos))
                        DDRB1trsmatch = DA_DRB1[SAAP_DRB1].nunique()
                            #for SAAP_DRB1 in DDRB1n:
                                #DDBR1trsmatch =DDRB1n.get(SAAP_DRB1)
                        if (DDRB1trsmatch == 1):
                            prob = 1
                            DAA_DRB1[subject_id][SAAP_DRB1]=prob
                        else:
                            prob = TRS_AA(valuesList_DDRB1)
                            DAA_DRB1[subject_id][SAAP_DRB1]=prob

                    DDQB1_sub =DGF_DQB1.get(subject_id)
                    valuesList_DDQB1 = list(DDQB1_sub.values())
                    for keys in DDQB1_sub:
                        DDQB1_list = keys.split('+')
                        DA_DQB1 = defaultdict(dict)
                        for DDQB1 in DDQB1_list:
                            for i in range(1,loci_range_value):
                                DDQB1loc = GENO_DQB1[0:4].split('*')[0]
                                DDQB1_ard_start = aa_mm.ard_start_pos[DDQB1loc]
                                DDQB1_ard_end = aa_mm.ard_end_pos[DDQB1loc]
                                aa_DDQB1_full = aa_mm.getAAsubstring(DDQB1,DDQB1_ard_start,DDQB1_ard_end)
                                aa_DDQB1= list(aa_DDQB1_full)
                                for pos in range(DQB1_ard_start,DQB1_ard_end):
                                    SAAP_DQB1 = str(DDQB1loc+"_"+str(pos))
                                    DADQB1 = aa_DDQB1[pos-1]
                                    DA_DQB1[DDQB1][SAAP_DQB1]=DADQB1
                    DA_DQB1 = pd.DataFrame.from_dict(DA_DQB1, orient = 'index')
                    #for column in DA_DQB1:
                        #DDQB1n = DA_DQB1.nunique(axis=0)
                        #DDQB1n = dict(DDQB1n)
                    for pos in range(DDQB1_ard_start,DDQB1_ard_end):
                        SAAP_DQB1 = str(DDQB1loc+"_"+str(pos))
                        DDQB1trsmatch = DA_DQB1[SAAP_DQB1].nunique()
                           #for SAAP_DQB1 in DDQB1n:
                                #DDQB1trsmatch =DDQB1n.get(SAAP_DQB1)
                        if (DDQB1trsmatch == 1):
                            prob = 1
                            DAA_DQB1[subject_id][SAAP_DQB1]=prob
                        else:
                            prob = TRS_AA(valuesList_DDQB1)
                            DAA_DQB1[subject_id][SAAP_DQB1]=prob
                else:
                    continue
            
    
  
            DAA_A = pd.DataFrame.from_dict(DAA_A, orient = 'index')
            DAA_A=DAA_A.mean(axis='index')
            Dtrs_aa = Dtrs_aa.append(pd.DataFrame([DAA_A], index=[pop] ))
            #print(Dtrs_aa)

            DAA_C = pd.DataFrame.from_dict(DAA_C, orient = 'index')
            #DC_dataframe = pd.DataFrame(DAA_C)
            #DC_design_matrix_filename = pathloc +"/AAgeno_DC_trs_5loc_aa_matrix_recip.csv"
            #DC_dataframe.to_csv(DC_design_matrix_filename,index= True,sep=",")
            DAA_C=DAA_C.mean(axis='index')
            Dtrs_aa = Dtrs_aa.append(pd.DataFrame([DAA_C], index=[pop] ))
            
            DAA_B = pd.DataFrame.from_dict(DAA_B, orient = 'index')
            DAA_B=DAA_B.mean(axis='index')
            Dtrs_aa = Dtrs_aa.append(pd.DataFrame([DAA_B], index=[pop] ))
            
            DAA_DRB1 = pd.DataFrame.from_dict(DAA_DRB1, orient = 'index')
            DAA_DRB1=DAA_DRB1.mean(axis='index')
            Dtrs_aa = Dtrs_aa.append(pd.DataFrame([DAA_DRB1], index=[pop] ))
            
            
            DAA_DQB1 = pd.DataFrame.from_dict(DAA_DQB1, orient = 'index')
            DAA_DQB1=DAA_DQB1.mean(axis='index')
            Dtrs_aa = Dtrs_aa.append(pd.DataFrame([DAA_DQB1], index=[pop] ))
            
            
            Dtrs_dataframe = pd.DataFrame(Dtrs_aa)
            Ddesign_matrix_filename = pathloc +"/trs_5loc_aa_matrix_donor.csv"
            Dtrs_dataframe.to_csv(Ddesign_matrix_filename,index= True,sep=",")
#debugging even distibution

exit()

