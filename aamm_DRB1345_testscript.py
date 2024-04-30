from aa_matching_msf_genie_MOD import *

def count_AA_Mismatches_Allele_DR(self, allele1_donor,allele2_donor,allele3_donor,
                                allele4_donor,allele1_recip,allele2_recip,
                                allele3_recip,allele4_recip,position):
        donor_homoz = 0
        if (allele1_donor == allele2_donor):
            donor_homoz+=1
        if (allele3_donor == allele4_donor):
            donor_homoz+=1
        return(donor_homoz)

        aa1_donor = self.getAAposition(allele1_donor,position)
        aa2_donor = self.getAAposition(allele2_donor,position)
        aa3_donor = self.getAAposition(allele3_donor,position)
        aa4_donor = self.getAAposition(allele4_donor,position)
        aa1_recip = self.getAAposition(allele1_recip,position)
        aa2_recip = self.getAAposition(allele2_recip,position)
        aa3_recip = self.getAAposition(allele3_recip,position)
        aa4_recip = self.getAAposition(allele4_recip,position)

        print(aa1_donor)
        print(aa2_donor)
        print(aa3_donor)
        print(aa4_donor)
        print(aa1_recip)
        print(aa2_recip)
        print(aa3_recip)
        print(aa4_recip)

        mm_count = self.count_AA_Mismatches_DR(aa1_donor,aa2_donor,aa3_donor,aa4_donor,aa1_recip,aa2_recip,aa3_recip,aa4_recip)

        if (donor_homoz == 1):
            mm_count-=1
        if (donor_homoz == 2):
            mm_count-=2

        return mm_count

allele1_donor = "DRB1*07:01"
allele2_donor = "DRB1*15:01"
allele3_donor = "DRB4*01:03"
allele4_donor = "DRB5*01:01"
allele1_recip = "DRB1*11:01"
allele2_recip = "DRB1*15:03"
allele3_recip = "DRB3*02:02"
allele4_recip = "DRB4*01:01"

aa_mm = AAMatch(dbversion=3420)
aam = aa_mm.count_AA_Mismatches_Allele_DR(allele1_donor,allele2_donor,allele3_donor,allele4_donor,allele1_recip,allele2_recip,allele3_recip,allele4_recip,179)

print(aam)
