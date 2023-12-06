#!/usr/bin/env python
from collections import defaultdict
import sys

# split pull and glid files into separate pops

# open pop combine file

pops_rollup = ['AFA','ASN','CAU','HIS','NAM','MLT']


# rollup US
for pop in pops_rollup:

    sys.stderr.write(pop + "\n")

    locilist = ['A','C','B','DRBX','DRB1','DQA1','DQB1','DPA1','DPB1']
    antigens = []
    for loc in locilist:
        greedy_antigens_filename = "./greedy/antigens." + pop + "." + loc + ".txt"
        greedy_antigens_file = open (greedy_antigens_filename, "r")
        for line in greedy_antigens_file:
            line = line.rstrip('\r\n')            
            antigens.append(line)

    # print (antigens)

    pull_rollup_filename = "./pull/pull.srtr." + pop + ".txt" 
    pull_rollup_file = open (pull_rollup_filename, "w")

    glid_rollup_filename = "./pull/glid.srtr." + pop + ".txt" 
    glid_rollup_file = open (glid_rollup_filename, "w")

    pull_filename = "./pull/pull.srtr.txt"
    pull_file = open(pull_filename, "r")

    glid_filename = "./pull/glid.srtr.txt"
    glid_file = open(glid_filename, "r")

    glid_greedy_filename = "./greedy/glid." + pop + ".txt"
    glid_greedy_file = open(glid_greedy_filename, "r")    

    glid_missing = {}

    for glid_line in glid_greedy_file:
        (glid,glstring) = glid_line.split(',')
        if (glid in ["0","999999999","888888888","55555555","11111111111","777777777"]):   
            glid_rollup_file.write(glid_line)
    
    glid_greedy_file.close()

    glids = {}
    for line in pull_file:
        line = line.rstrip('\r\n')
        (id,rollup_race,glid_a,glid_b,glid_c,glid_drb1,glid_drbx,glid_dqa1,glid_dqb1,glid_dpa1,glid_dpb1) = line.split(',')
        # print (id)
        if (rollup_race != pop):
            continue

        # rename glid 0 for most loci
        if (glid_dqb1 == "0"):
            glid_dqb1 = "999999999"
        if (glid_drbx == "0"):
            glid_drbx = "888888888"
        if (glid_dqa1 == "0"):
            glid_dqa1 = "55555555"
        if (glid_dpa1 == "0"):
            glid_dpa1 = "11111111111"
        if (glid_dpb1 == "0"):
            glid_dpb1 = "777777777"

        glids[glid_a] = 1
        glids[glid_b] = 1
        glids[glid_c] = 1
        glids[glid_drb1] = 1
        glids[glid_drbx] = 1
        glids[glid_dqa1] = 1
        glids[glid_dqb1] = 1
        glids[glid_dpa1] = 1
        glids[glid_dpb1] = 1

        glidlist = ":".join([glid_a,glid_c,glid_b,glid_drbx,glid_drb1,glid_dqa1,glid_dqb1,glid_dpa1,glid_dpb1])
        pull_rollup_file.write(id + "," + glidlist + "\n")

    # print (glids)
    
    for glid_line in glid_file:
        # print (glid_line)
        glid_line = glid_line.rstrip('\r\n')
        (glid,glstring) = glid_line.split(',')

        # skip glids not in cases
        if (glid not in glids):
            continue

        # reformat glstring to "+" and "|" only
        
        genos_full = []

        # make intermediate genotype list
        genos_inter = glstring.split('|')

        for geno in genos_inter:
            loci = geno.split('+')
            # print (loci)
            if (len(loci) == 1):
                loci.append(loci[0])
            if (loci[0] == ""):
                loci[0] = loci[1] # 73354:+DRB1*15:01
            # strip out DRB345
            l1_dr = loci[0].split('~')
            l2_dr = loci[1].split('~')
            l1_array = l1_dr[0].split('/')
            l2_array = l2_dr[0].split('/')

            for loctyp1 in l1_array:
                for loctyp2 in l2_array:
                    if ((loctyp1 in antigens) & (loctyp2 in antigens)):
                        genos_full.append(loctyp1 + "+" + loctyp2)

        if ((len(genos_full) == 0)):

            # sys.stderr.write("Don't reduce GLID: " + glid + "," + glstring + "\n")

            # should be short with rare alleles
            for geno in genos_inter:
                loci = geno.split('+')
                # print (loci)
                if (len(loci) == 1):
                    loci.append(loci[0])
                if (loci[0] == ""):
                    loci[0] = loci[1] # 73354:+DRB1*15:01
                # strip out DRB345
                l1_dr = loci[0].split('~')
                l2_dr = loci[1].split('~')
                l1_array = l1_dr[0].split('/')
                l2_array = l2_dr[0].split('/')

                for loctyp1 in l1_array:
                    for loctyp2 in l2_array:
                        genos_full.append(loctyp1 + "+" + loctyp2)


        glstring = "|".join(genos_full)

        # Update genotype list to consider missing DRB345 or DRB4 (due to lack of intent)
        if ((glstring.startswith("DRB3")) | (glstring.startswith("DRB4")) | (glstring.startswith("DRB5")) | (glstring.startswith("DRBX"))): 
            homo_DRB345 = 0   # 1 if homozygous typing found
            DRB345_typs = {} # dictionary of DRBX typings
            for geno in genos_full:
                # sys.stderr.write(geno + "\n")
                (loctyp1,loctyp2) = geno.split("+")
                if (loctyp1 == loctyp2):
                    homo_DRB345 = 1
                DRB345_typs[loctyp1] = 1
            if (homo_DRB345 == 1):
                for DRB345_loctyp1 in DRB345_typs:
                    geno = DRB345_loctyp1 + "+DRBX*NNNN"
                    genos_full.append(geno)
                    (loc1, typ1) = DRB345_loctyp1.split("*")

                    if (loc1 == "DRB3"):
                        geno = DRB345_loctyp1 + "+DRB4*01:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB5*01:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB5*01:02"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB5*02:02"
                        genos_full.append(geno)
                    
                    if (loc1 == "DRB4"):
                        geno = DRB345_loctyp1 + "+DRB3*01:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB3*02:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB3*02:02"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB3*03:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB5*01:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB5*01:02"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB5*02:02"
                        genos_full.append(geno)

                    if (loc1 == "DRB5"):
                        geno = DRB345_loctyp1 + "+DRB3*01:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB3*02:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB3*02:02"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB3*03:01"
                        genos_full.append(geno)
                        geno = DRB345_loctyp1 + "+DRB4*01:01"
                        genos_full.append(geno)

            # sys.stderr.write("DRB345 Pre:" + glstring + "\n")
            glstring = "|".join(genos_full)
            # sys.stderr.write("DRB345 Post:" + glstring + "\n")

        if glid in glids:
            # print (glid)
            glid_rollup_file.write(glid + "," + glstring + "\n")


exit(0)  