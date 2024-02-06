#!/usr/bin/env python

import gzip
from collections import defaultdict
import hlagenie
genie = hlagenie.init("3510")

loci = ["A", "B", "C", "DRB345", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

full_start_pos = {
    "A" : 1,
    "B" : 1,
    "C" : 1,
    "DRB345" : 1,
    "DRB1" : 1,
    "DQA1" : 1,
    "DQB1" : 1,
    "DPA1" : 1,
    "DPB1" : 1,
}

full_end_pos = {
    "A": 341,
    "B": 338,
    "C": 342,
    "DRB1": 237,
    "DRB345": 237,
    "DQA1": 232,
    "DQB1": 229,  # increased by 1
    "DPA1": 229,
    "DPB1": 229,
}


# Utility function to create dictionary
def multi_dict(K, type):
    if K == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: multi_dict(K-1, type))

pops = ['AFA', 'ASN', 'CAU', 'HIS', 'MLT', 'NAM']
happair_probs = {}  # HLA probability distribution
happair_hla = {}  # HLA haplotype pair distribution
happair_id_total = {}  # cumulative genotype frequency total
subject_ID_ethnicity_study = {}
for pop in pops:
    # load SRTR HapLogic imputation output file and load probabilities for all subjects
    impute_outfilename = "./impute.srtr." + pop + ".csv.gz"
    impute_outfile = gzip.open(impute_outfilename, "rt")
    # compute cumulative genotype frequency totals per subject

    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, freq) = line.split(',')
        if (subject_id == "PX_ID"):  # skip header row
            continue

        happair_freq = 0

        subject_ID_ethnicity_study[subject_id] = pop

        if (hap1 == hap2):
            happair_freq = float(freq)
        else:
            happair_freq = 2 * float(freq)

        if subject_id not in happair_id_total:
            happair_id_total[subject_id] = 0
        happair_id_total[subject_id] += float(happair_freq)

        # print (subject_id + " " + hap1 + "+" + hap2 + " " + str(happair_id_total[subject_id]))
    impute_outfile.close()

    impute_outfile = gzip.open(impute_outfilename, "rt")

    # compute probabilties for each haplotype pair
    for line in impute_outfile:
        (subject_id, rank, hap1, hap2, freq) = line.split(',')
        if (subject_id == "PX_ID"):  # skip header row
            continue

        happair = hap1 + "+" + hap2

        happair_freq = 0
        happair_total = 0
        happair_freq = float(freq)

        if subject_id not in happair_probs:
            happair_probs[subject_id] = list()
        happair_probs[subject_id].append(happair_freq / happair_id_total[subject_id])
        if subject_id not in happair_hla:
            happair_hla[subject_id] = list()
        happair_hla[subject_id].append(happair)
        # print (subject_id + " " + happair + " " + str(happair_freq))


# compute probabilities for each AA-level genotype
TRS_dict = multi_dict(3, float)  # TRS per subject per locus per position
subject_counter = 0
for subject_id in happair_hla:

    # starttime = timeit.default_timer()

    happair_index = 0
    AA_geno_probs = multi_dict(3, float)  # AA-level genotype list and probs per locus position
    for happair in happair_hla[subject_id]:
        # print (subject_id)
        # print ("Haplotype Pair: " + happair)
        (hap1, hap2) = happair.split('+')
        # print(hap1)
        (a1, c1, b1, dr345_1, drb1_1, dqa1_1, dqb1_1, dpa1_1, dpb1_1) = hap1.split('~')
        (a2, c2, b2, dr345_2, drb1_2, dqa1_2, dqb1_2, dpa1_2, dpb1_2) = hap2.split('~')

        if a1 == "A*11:52":
            a1 = "A*11:52Q"
        if a2 == "A*11:52":
            a2 = "A*11:52Q"
        if a1 == "A*23:19Q":
            a1 = "A*23:19N"
        if a2 == "A*23:19Q":
            a2 = "A*23:19N"

        if b1 == "B*13:08Q":
            b1 = "B*13:08"
        if b2 == "B*13:08Q":
            b2 = "B*13:08"
        if b1 == "B*15:22":
            b1 = "B*35:43"
        if b2 == "B*15:22":
            b2 = "B*35:43"

        if c1 == "C*15:20":
            c1 = "C*15:27"
        if c2 == "C*15:20":
            c2 = "C*15:27"
        if c1 == "C*03:12":
            c1 = "C*03:19"
        if c2 == "C*03:12":
            c2 = "C*03:19"
        if c1 == "C*03:23":
            c1 = "C*03:23N"
        if c2 == "C*03:23":
            c2 = "C*03:23N"
        if c1 == "C*05:02":
            c1 = "C*05:09"
        if c2 == "C*05:02":
            c2 = "C*05:09"
        if c1 == "C*13:01":
            c1 = "C*12:02"
        if c2 == "C*13:01":
            c2 = "C*12:02"

        if drb1_1 == "DRB1*07:02":
            drb1_1 = "DRB1*07:01"
        if drb1_2 == "DRB1*07:02":
            drb1_2 = "DRB1*07:01"

        if dqa1_1 == "DQA1*01:07":
            dqa1_1 = "DQA1*01:07Q"
        if dqa1_2 == "DQA1*01:07":
            dqa1_2 = "DQA1*01:07Q"

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
            if (locus == "DRB345"):
                allele1 = dr345_1
                allele2 = dr345_2
            if (locus == "DRB1"):
                allele1 = drb1_1
                allele2 = drb1_2
            if (locus == "DQA1"):
                allele1 = dqa1_1
                allele2 = dqa1_2
            if (locus == "DQB1"):
                allele1 = dqb1_1
                allele2 = dqb1_2
            if (locus == "DPA1"):
                allele1 = dpa1_1
                allele2 = dpa1_2
            if (locus == "DPB1"):
                allele1 = dpb1_1
                allele2 = dpb1_2

            if subject_id not in happair_probs:
                happair_probs[subject_id] = list()

            # get appropriate position range per locus
            # lots of incomplete sequences in IMGT/HLA that are filled in - use only ARD positions
            # for position in DRB1_AA_positions_selected:
            for position in range(full_start_pos[locus], full_end_pos[locus]):
                # print (locus)
                # print ("DRB1 allele-level genotypes: " + drb1_1 + "+" + drb1_2)
                if allele1 == "DRBX*NNNN":
                    AA1 = 'None'
                else:
                    AA1 = genie.getAA(allele1, position)

                if allele2 == "DRBX*NNNN":
                    AA2 = 'None'
                else:
                    AA2 = genie.getAA(allele2, position)

                (AA1, AA2) = sorted([AA1, AA2])  # alpha sort positions
                AA_geno = AA1 + "+" + AA2
                # print ("AA position: " + str(position))
                # print ("AA-level genotype: " + AA_geno)
                # print ("Hap Pair Prob: " + str(happair_prob))
                AA_geno_probs[locus][position][AA_geno] += happair_prob

    for locus in loci:
        # for locus in loci_selected:
        # print (locus)
        # for position in DRB1_AA_positions_selected:
        for position in range(full_start_pos[locus], full_end_pos[locus]):
            # print ("AA position: " + str(position))
            TRS = 0
            for AA_geno in AA_geno_probs[locus][position]:
                AA_geno_prob = AA_geno_probs[locus][position][AA_geno]
                AA_geno_prob = round(AA_geno_prob, 10)
                # print ("Amino acid level genotype and probs: " + AA_geno + ": " + str(AA_geno_prob))
                TRS = TRS + (AA_geno_prob * AA_geno_prob)
                TRS_increment = AA_geno_prob * AA_geno_prob
                # print ("AA genotype contribution to TRS: " + str(TRS_increment))
            TRS_dict[subject_id][locus][position] = TRS

    # print("The time difference is :", timeit.default_timer() - starttime)
    subject_counter += 1
    if ((subject_counter % 10000) == 0):
        print("Number of subjects with TRS computed: " + str(subject_counter))

# print out summary of all TRS values
TRS_average = multi_dict(4, float)  # Average TRS per donor/recip per race per locus per position



nsubjects = len(subject_ID_ethnicity_study)
print ("Number of subjects in the impute.srtr.*.csv.gz file: " + str(nsubjects))

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
print ("Number of MLT donors: " + str(nsubject_donor_ethnicity["MLT"]))
print ("Number of MLT recips: " + str(nsubject_recip_ethnicity["MLT"]))

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
        for position in range(full_start_pos[locus], full_end_pos[locus]):
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
TRS_average_out_file = open(TRS_average_out_filename, "w")
TRS_average_out_file.write("Subject_Type,Ethnicity,Locus,AA_Position,TRS_Average\n")
for donor_or_recip in ["SUBJECT", "DONOR", "RECIP"]:
    for ethnicity in ethnicity_list:
        # for locus in loci_selected:
        for locus in loci:
            for position in range(full_start_pos[locus], full_end_pos[locus]):
                # for position in DRB1_AA_positions_selected:
                print(donor_or_recip + " " + ethnicity + " " + locus + " " + str(position))

                # UNK recip has 0 cases, so TRS should be NA
                TRS = ""
                if (not TRS_average[donor_or_recip][ethnicity][locus][position]):
                    TRS = "NA"
                else:
                    TRS = str(round(TRS_average[donor_or_recip][ethnicity][locus][position], 8))

                print("Average TRS: " + TRS)
                TRS_out_string = ",".join([donor_or_recip, ethnicity, locus, str(position), TRS])
                TRS_average_out_file.write(TRS_out_string + "\n")

# TODO - create separate table per locus
# rows in each tables are donor/recip-ethnicity combo
# columns are positions


exit()
