from aa_matching_msf_genie_MOD import *

# test DR amino acid mismatch computation, considering both DRB1 and DRB3/4/5 loci as homologous

allele1_donor = "DRB1*07:01"
allele2_donor = "DRB1*15:01"
allele3_donor = "DRB4*01:03"
allele4_donor = "DRB5*01:01"
allele1_recip = "DRB1*11:01"
allele2_recip = "DRB1*15:03"
allele3_recip = "DRB3*02:02"
allele4_recip = "DRB4*01:01"

position = 13 # most variable position in DRB1

aa_mm = AAMatch(dbversion=3420)
aam = aa_mm.count_AA_Mismatches_Allele_DR(allele1_donor,allele2_donor,allele3_donor,allele4_donor,allele1_recip,allele2_recip,allele3_recip,allele4_recip,position)

print(aam)
