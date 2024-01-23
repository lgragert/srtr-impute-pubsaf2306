#!/usr/bin/env python

from pathlib import Path
import pandas as pd
from collections import defaultdict
import gzip

# Compute typing resolution score per population and per position
pathloc = str(Path('.'))
population_set = "ALL"


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


pops_NEMO = ['AFA','API','CAU','HIS','NAM','AAFA','AFB','CARB','SCAMB','AINDI','FILII','HAWI','JAPI','KORI','NCHI','SCSEAI','VIET','EURCAU','MENAFC','CARHIS','MSWHIS','SCAHIS','AISC','ALANAM','AMIND','CARIBI']
pops_US = ['AFA','API','CAU','HIS','NAM','DEC','MLT','OTH','UNK','AAFA','AFB','CARB','SCAMB','AINDI','FILII','HAWI','JAPI','KORI','NCHI','SCSEAI','VIET','EURCAU','MENAFC','CARHIS','MSWHIS','SCAHIS','AISC','ALANAM','AMIND','CARIBI']
pops_GLOBAL = ['US','DE','BR','IL','GB','PL','CN','TR','IN','ES','JP','CA','PT','NL','IT','FR','TH','AR','GR','SE','CY','CH','AU','SG','CL','CZ','SA','MX','AT','BE','HR','DK','RO','FI','IR','ZA','AM','NO','RU','SI','IE','SK','LT','NZ','RS','HU','BG','MK','UY','NG','GLOBAL']
pops_CPRA = ['AFA','CAU','HIS','NAM','ASN','HPI','MLT','HISETH']
pops_ALL = ['AFA','ASN','CAU','HIS','NAM']  # ['NAM']
patient = ['recip','donor']

trs = []

if population_set == "ALL":
    pops = pops_ALL

print("Pop,TRS_A,TRS_C,TRS_B,TRS_DRB1,TRS_DRB345,TRS_DQA1,TRS_DQB1,TRS_DPA1,TRS_DPB1")


for pop in pops:

    print(pop)
    impute_outfilename = pathloc + "/impute.srtr." + pop + ".csv.gz"
    print(impute_outfilename)
    impute_outfile = gzip.open(impute_outfilename, "rt")

    # Overall frequency already calculated for haplo pairs in 9loc impute data

    # compute cumulative genotype frequency totals per subject for later renormalization
    happair_id_total = {}  # cumulative genotype frequency total

    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, freq) = line.split(',')

        happair_freq = 0
        happair_freq = float(freq)

        if subject_id not in happair_id_total:
            happair_id_total[subject_id] = 0
        happair_id_total[subject_id] += float(happair_freq)

    impute_outfile.close()

    ndonor = len(happair_id_total)

    # reopen impute file for 2nd pass
    impute_outfile = gzip.open(impute_outfilename, "rt")

    # compute probabilties for each haplotype pair
    happair_probs = {}  # HLA probability distribution
    happair_hla = {}  # HLA haplotype pair distribution

    # Dicts to store allele level probs for genotypes
    GF_A = defaultdict(dict)
    GF_B = defaultdict(dict)
    GF_C = defaultdict(dict)
    GF_DRB345 = defaultdict(dict)
    GF_DRB1 = defaultdict(dict)
    GF_DQA1 = defaultdict(dict)
    GF_DQB1 = defaultdict(dict)
    GF_DPA1 = defaultdict(dict)
    GF_DPB1 = defaultdict(dict)

    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, freq) = line.split(',')

        (A_1, C_1, B_1, DRB345_1, DRB1_1, DQA1_1, DQB1_1, DPA1_1, DPB1_1) = hap1.split('~')
        (A_2, C_2, B_2, DRB345_2, DRB1_2, DQA1_2, DQB1_2, DPA1_2, DPB1_2) = hap2.split('~')

        happair = hap1 + "+" + hap2
        happair_freq = float(freq)

        # sort strings
        (A_1, A_2) = sorted([A_1, A_2])
        (B_1, B_2) = sorted([B_1, B_2])
        (C_1, C_2) = sorted([C_1, C_2])
        (DRB345_1, DRB345_2) = sorted([DRB345_1, DRB345_2])
        (DRB1_1, DRB1_2) = sorted([DRB1_1, DRB1_2])
        (DQA1_1, DQA1_2) = sorted([DQA1_1, DQA1_2])
        (DQB1_1, DQB1_2) = sorted([DQB1_1, DQB1_2])
        (DPA1_1, DPA1_2) = sorted([DPA1_1, DPA1_2])
        (DPB1_1, DPB1_2) = sorted([DPB1_1, DPB1_2])

        GENO_A = '+'.join([A_1, A_2])
        GENO_C = '+'.join([C_1, C_2])
        GENO_B = '+'.join([B_1, B_2])
        GENO_DRB345 = '+'.join([DRB345_1, DRB345_2])
        GENO_DRB1 = '+'.join([DRB1_1, DRB1_2])
        GENO_DQA1 = '+'.join([DQA1_1, DQA1_2])
        GENO_DQB1 = '+'.join([DQB1_1, DQB1_2])
        GENO_DPA1 = '+'.join([DPA1_1, DPA1_2])
        GENO_DPB1 = '+'.join([DPB1_1, DPB1_2])

        prob = happair_freq / happair_id_total[subject_id]

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

        if GENO_DRB345 not in GF_DRB345[subject_id]:
            GF_DRB345[subject_id][GENO_DRB345] = prob
        else:
            GF_DRB345[subject_id][GENO_DRB345] = GF_DRB345[subject_id][GENO_DRB345] + prob

        if GENO_DRB1 not in GF_DRB1[subject_id]:
            GF_DRB1[subject_id][GENO_DRB1] = prob
        else:
            GF_DRB1[subject_id][GENO_DRB1] = GF_DRB1[subject_id][GENO_DRB1] + prob

        if GENO_DQA1 not in GF_DQA1[subject_id]:
            GF_DQA1[subject_id][GENO_DQA1] = prob
        else:
            GF_DQA1[subject_id][GENO_DQA1] = GF_DQA1[subject_id][GENO_DQA1] + prob

        if GENO_DQB1 not in GF_DQB1[subject_id]:
            GF_DQB1[subject_id][GENO_DQB1] = prob
        else:
            GF_DQB1[subject_id][GENO_DQB1] = GF_DQB1[subject_id][GENO_DQB1] + prob

        if GENO_DPA1 not in GF_DPA1[subject_id]:
            GF_DPA1[subject_id][GENO_DPA1] = prob
        else:
            GF_DPA1[subject_id][GENO_DPA1] = GF_DPA1[subject_id][GENO_DPA1] + prob

        if GENO_DPB1 not in GF_DPB1[subject_id]:
            GF_DPB1[subject_id][GENO_DPB1] = prob
        else:
            GF_DPB1[subject_id][GENO_DPB1] = GF_DPB1[subject_id][GENO_DPB1] + prob

    (TRS_A_subject, TRS_A_mean) = TRS_locus(GF_A)
    (TRS_B_subject, TRS_B_mean) = TRS_locus(GF_B)
    (TRS_C_subject, TRS_C_mean) = TRS_locus(GF_C)
    (TRS_DRB345_subject, TRS_DRB345_mean) = TRS_locus(GF_DRB345)
    (TRS_DRB1_subject, TRS_DRB1_mean) = TRS_locus(GF_DRB1)
    (TRS_DQA1_subject, TRS_DQA1_mean) = TRS_locus(GF_DQA1)
    (TRS_DQB1_subject, TRS_DQB1_mean) = TRS_locus(GF_DQB1)
    (TRS_DPA1_subject, TRS_DPA1_mean) = TRS_locus(GF_DPA1)
    (TRS_DPB1_subject, TRS_DPB1_mean) = TRS_locus(GF_DPB1)

    trs.append({
        'Population': pop,
        'TRS_A_mean': TRS_A_mean,
        'TRS_C_mean': TRS_C_mean,
        'TRS_B_mean': TRS_B_mean,
        'TRS_DRB1_mean': TRS_DRB1_mean,
        'TRS_DRW1_mean': TRS_DRB345_mean,
        'TRS_DQA1_mean': TRS_DQA1_mean,
        'TRS_DQB1_mean': TRS_DQB1_mean,
        'TRS_DPA1_mean': TRS_DPA1_mean,
        'TRS_DPB1_mean': TRS_DPB1_mean
    })

    trs_dataframe = pd.DataFrame(trs)
    design_matrix_filename = pathloc + "/trs_9loc_matrix.csv"
    trs_dataframe.to_csv(design_matrix_filename, index=False, sep=",")
    impute_outfile.close()

exit()
