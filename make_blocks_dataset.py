#!/usr/bin/env python
from collections import defaultdict
import sys
import os
import time
import gzip

# Make Blocks for imputation pipeline
# Merges imputation output with pull to make two-locus blocks

dataset = sys.argv[1]
pop = sys.argv[2]
loc1 = sys.argv[3]
loc2 = sys.argv[4]
pull_folder = sys.argv[5]
input_freqs_global = sys.argv[6] # 1 - GLOBAL subdirectory, 0 - NON-GLOBAL
loc_folder_out = sys.argv[7]
impute_threshold = sys.argv[8]  # default 0.01

start_time = time.time()

sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)

sys.stderr.write("Make Blocks: " + loc1 + " " + loc2 + "\n")

# locus positions in GLID array in pull file
loc_position = {}
loc_position["A"] = 0
loc_position["C"] = 1
loc_position["B"] = 2
loc_position["DRBX"] = 3
# loc_position["DRB3"] = 3
# loc_position["DRB4"] = 3
# loc_position["DRB5"] = 3
loc_position["DRB1"] = 4
loc_position["DQA1"] = 5
loc_position["DQB1"] = 6
loc_position["DPA1"] = 7
loc_position["DPB1"] = 8

# get GLIDs from Cases GLID file
glstring = {} # stores GL String for each GLID
glid_filename = "./" + pull_folder + "/glid." + dataset + "." + pop + ".txt"
glid_file = open(glid_filename,'r')
for line in glid_file:
    line = line.rstrip('\r\n')
    (glid,gl) = line.split(',')
    (loc,*string) = gl.split('*')

    # handle DRB345
    if (loc == "DRB3") | (loc == "DRB4") | (loc == "DRB5"):
        loc = "DRBX"

    # won't find loci with "~"
    if (loc == loc1) | (loc == loc2):
        glstring[glid] = gl
glid_file.close()


# map tilde to filename
loc_combo = {}
if (input_freqs_global == "0"):
    loc_combo["C~B"] = "CB"
    loc_combo["A~C~B"] = "ACB"
    loc_combo["DRBX~DRB1"] = "DRBXDRB1"
    loc_combo["DRBX~DRB1~DQB1"] = "DRBXDRB1DQB1"
    loc_combo["DRBX~DRB1~DQA1~DQB1"] = "DRBXDRB1DQA1DQB1"
    loc_combo["DRBX~DRB1~DQA1~DQB1~DPB1"] = "DRBXDRB1DQA1DQB1DPB1"
    loc_combo["DRBX~DRB1~DQA1~DQB1~DPA1~DPB1"] = "DRBXDRB1DQA1DQB1DPA1DPB1"

    loc_combo["DRBX~DRB1~DQB1~DPB1"] = "DRBXDRB1DQB1DPB1"
    loc_combo["DQA1~DQB1"] = "DQA1DQB1"
    loc_combo["DPA1~DPB1"] = "DPA1DPB1"

    loc_combo["DRB1~DQB1"] = "DRB1DQB1"
    loc_combo["DRB1~DQA1~DQB1"] = "DRB1DQA1DQB1"
else:
    loc_combo["C~B"] = "CB_GLOBAL"
    loc_combo["A~C~B"] = "ACB_GLOBAL"
    loc_combo["DRBX~DRB1"] = "DRBXDRB1_GLOBAL"
    loc_combo["DRBX~DRB1~DQB1"] = "DRBXDRB1DQB1_GLOBAL"
    loc_combo["DRBX~DRB1~DQA1~DQB1"] = "DRBXDRB1DQA1DQB1_GLOBAL"
    loc_combo["DRBX~DRB1~DQA1~DQB1~DPB1"] = "DRBXDRB1DQA1DQB1DPB1_GLOBAL"
    loc_combo["DRBX~DRB1~DQA1~DQB1~DPA1~DPB1"] = "DRBXDRB1DQA1DQB1DPA1DPB1_GLOBAL"

    loc_combo["DRBX~DRB1~DQB1~DPB1"] = "DRBXDRB1DQB1DPB1_GLOBAL"
    loc_combo["DQA1~DQB1"] = "DQA1DQB1_GLOBAL"
    loc_combo["DPA1~DPB1"] = "DPA1DPB1_GLOBAL"

    loc_combo["DRB1~DQB1"] = "DRB1DQB1_GLOBAL"
    loc_combo["DRB1~DQA1~DQB1"] = "DRB1DQA1DQB1_GLOBAL"

# if locus name has "~", use imputation output file
if ("~" in loc1):
    loc_folder = loc_combo[loc1]
    impute_filename = "./" + loc_folder + "/impute." + dataset + "." + pop + ".csv.gz"
    impute_file = gzip.open(impute_filename,'rt')
    freqnorm_highest = {}
    imputed_genos = {}
    for line in impute_file:
        line = line.rstrip('\r\n')
        (id, rank, h1, h2, freqnorm) = line.split(',')

        # skip less probable genotypes
        if (id not in freqnorm_highest): # set freq of zero for ID
            freqnorm_highest[id] = 0
        if (float(freqnorm) > float(freqnorm_highest[id])):
            freqnorm_highest[id] = freqnorm
        if ((float(freqnorm) < float(impute_threshold)) & (float(freqnorm_highest[id]) > float(impute_threshold))):
            continue

        # sort genotypes
        if (h1 > h2):
            (h2,h1) = (h1,h2)
        
        geno = "+".join([h1,h2])
    
        if (id not in imputed_genos):
            imputed_genos[id] = {}
        imputed_genos[id][geno] = 1

        # print (id + " " + geno)
    
    impute_file.close()

    GLstrings_block = {} # GLstrings made from block
    GLID_block1 = {} # GLID make from block for each donor
    GLID_counter = 1000000 # starting number for new GLIDs
    # create new GLIDs for genotype lists
    for id in imputed_genos:

      # sort keys to limit unique GLIDs created
      GL = []
      for geno in imputed_genos[id]:
        GL.append(geno)

      gl = "|".join(GL)

      # if GLstring is new, create new GLID
      if (gl not in GLstrings_block):
        glstring[GLID_counter] = gl
        GLID_block1[id] = GLID_counter
        GLstrings_block[gl] = GLID_counter
        GLID_counter += 1
      else:
        GLID_block1[id] = GLstrings_block[gl]

      # print (id + GLID_block[id] + '\n')


# if locus name has "~", use imputation output file
if ("~" in loc2):
    loc_folder = loc_combo[loc2]
    impute_filename = "./" + loc_folder + "/impute." + dataset + "." + pop + ".csv.gz"
    impute_file = gzip.open(impute_filename,'rt')
    freqnorm_highest = {}
    imputed_genos = {}
    for line in impute_file:
        line = line.rstrip('\r\n')
        (id, rank, h1, h2, freqnorm) = line.split(',')

        # skip less probable genotypes
        if (id not in freqnorm_highest): # set freq of zero for ID
            freqnorm_highest[id] = 0
        if (float(freqnorm) > float(freqnorm_highest[id])):
            freqnorm_highest[id] = freqnorm
        # print (id + "," + freqnorm + "," + impute_threshold + "," + freqnorm_highest[id])
        if ((float(freqnorm) < float(impute_threshold)) & (float(freqnorm_highest[id]) > float(impute_threshold))):
            continue

        # sort genotypes
        if (h1 > h2):
            (h2,h1) = (h1,h2)
        
        geno = "+".join([h1,h2])
    
        if (id not in imputed_genos):
            imputed_genos[id] = {}
        imputed_genos[id][geno] = 1
        
        # print (id + " " + geno)
        # exit()
    
    impute_file.close()

    GLstrings_block = {} # GLstrings made from block
    GLID_block2 = {} # GLID make from block for each donor
    GLID_counter = 2000000 # starting number for new GLIDs
    # create new GLIDs for genotype lists
    for id in imputed_genos:

      # sort keys to limit unique GLIDs created
      GL = []
      for geno in imputed_genos[id]:
        GL.append(geno)

      gl = "|".join(GL)

      # if GLstring is new, create new GLID
      if (gl not in GLstrings_block):
        glstring[GLID_counter] = gl
        GLID_block2[id] = GLID_counter
        GLstrings_block[gl] = GLID_counter
        GLID_counter += 1
      else:
        GLID_block2[id] = GLstrings_block[gl]

      # print (str(id) + " " + str(GLID_block2[id]))

# Load Pull file
pull_filename = "./" + pull_folder + "/pull." + dataset + "." + pop + ".txt"

ID = {} # stores GLIDs per ID
GLIDlist_count = defaultdict(int) # stores count for set of GLIDs
pull_file = open(pull_filename,'r')
for line in pull_file:
    line = line.rstrip('\r\n')
    (id,glidlist) = line.split(',')
    glids = glidlist.split(':')

    new_glids = []

    if (loc1 in loc_position):
        loc1_position = loc_position[loc1]
        new_glids.append(glids[loc1_position])
    else:     # add GLIDs from imputation output
        new_glids.append(str(GLID_block1[id]))

    if (loc2 in loc_position):
        loc2_position = loc_position[loc2]
        new_glids.append(glids[loc2_position])
    else:     # add GLIDs from imputation output
        new_glids.append(str(GLID_block2[id]))

    glidlist = ":".join(new_glids)
    ID[id] = glidlist
    GLIDlist_count[glidlist] += 1

pull_file.close()


#  load DRBX~DRB1 association rules file
if ((loc1 == "DRBX") & (loc2 == "DRB1")):

    DRBX_DRB1_assoc = {}
    DRBX_DRB1_assoc_filename = "./cfg/DRBX_DRB1_assoc.cfg"
    DRBX_DRB1_assoc_file = open (DRBX_DRB1_assoc_filename,'r')
    for line in DRBX_DRB1_assoc_file:
        line = line.rstrip('\r\n')
        (drb1,drbx) = line.split(',')
        if (drb1 == "DRB1"):
            continue
        DRBX_DRB1_assoc[drb1] = drbx
        # print STDERR "$drb1 $drbx $DRBX_DRB1_assoc{$drb1}\n";

    DRBX_DRB1_assoc_file.close()

    # make all possible combo of DRBX genos with DRB1 genos
    genos_id_DRBX = defaultdict(dict) # final list of DRBX genotypes
    genos_id_DRB1 = defaultdict(dict) # final list of DRB1 genotypes
    restricted_DRBX_DRB1_filename = "./DRBXDRB1/restricted.DRBXDRB1." + dataset + "." + pop + ".txt"
    restricted_DRBX_DRB1_file = open (restricted_DRBX_DRB1_filename,'w')

    for id in ID:
        glidlist = ID[id]
        (glid_DRBX,glid_DRB1) = glidlist.split(":")
        glstring_DRBX = glstring[glid_DRBX]
        glstring_DRB1 = glstring[glid_DRB1]
        # print "DRBX: $glstring_DRBX DRB1: $glstring_DRB1\n";
        genos_DRBX = glstring_DRBX.split("|")
        genos_DRB1 = glstring_DRB1.split("|")
        genos_DRBX_restricted = [] # genos with DRBX assoc rules enforced
        genos_DRB1_restricted = [] # genos with DRBX assoc rules enforced
        genos_DRBX_all = [] # all genotypes - fallback if all restricted
        genos_DRB1_all = [] # all genotypes - fallback if all restricted
        genos_DRBXDRB1_all = []
        # genos_id_DRBX[id] = 1
        # genos_id_DRB1[id] = 1
        for geno_DRBX in genos_DRBX:

            (drbx_1,drbx_2) = geno_DRBX.split("+")
            if (drbx_1 > drbx_2):
                (drbx_1,drbx_2) = (drbx_2,drbx_1)
            geno_DRBX_sorted = "+".join([drbx_1,drbx_2])
            drbx_loc1 = drbx_1[0:4]
            drbx_loc2 = drbx_2[0:4]
            for geno_DRB1 in genos_DRB1:
                (drb1_1,drb1_2) = geno_DRB1.split("+")
                (drb1_loc1,drb1_typ1) = drb1_1.split("*")
                (drb1_loc2,drb1_typ2) = drb1_2.split("*")
                drb1_1_2dig = drb1_typ1[0:2]
                drb1_2_2dig = drb1_typ2[0:2]

                # phase 1
                geno_DRBXDRB1_phase1_1 = drbx_1 + "~" + drb1_1
                geno_DRBXDRB1_phase1_2 = drbx_2 + "~" + drb1_2
                
                if (geno_DRBXDRB1_phase1_1 > geno_DRBXDRB1_phase1_2):
                    (geno_DRBXDRB1_phase1_1,geno_DRBXDRB1_phase1_2) = (geno_DRBXDRB1_phase1_2,geno_DRBXDRB1_phase1_1)

                geno_phase1 = geno_DRBXDRB1_phase1_1 + "+" + geno_DRBXDRB1_phase1_2
                genos_DRBX_all.append(geno_DRBX_sorted)
                genos_DRB1_all.append(geno_DRB1)
                genos_DRBXDRB1_all.append(geno_phase1)

                # check for assoc rules
                # print ("Phase 1: DRB1 " + drb1_1_2dig + " DRBX " + drbx_loc1 + " DRB1 " + drb1_2_2dig +
                #         " DRBX " +  drbx_loc2)
                if ((DRBX_DRB1_assoc[drb1_1_2dig] == drbx_loc1) &
                    (DRBX_DRB1_assoc[drb1_2_2dig] == drbx_loc2)):
                    genos_DRBX_restricted.append(geno_DRBX_sorted)
                    genos_DRB1_restricted.append(geno_DRB1)
                    # print ("Accepted!")

                # phase 2
                geno_DRBXDRB1_phase2_1 = drbx_1 + "~" + drb1_2
                geno_DRBXDRB1_phase2_2 = drbx_2 + "~" + drb1_1

                if (geno_DRBXDRB1_phase2_1 > geno_DRBXDRB1_phase2_2):
                    (geno_DRBXDRB1_phase2_1,geno_DRBXDRB1_phase2_2) = (geno_DRBXDRB1_phase2_2,geno_DRBXDRB1_phase2_1)

                geno_phase2 = geno_DRBXDRB1_phase2_1 + "+" + geno_DRBXDRB1_phase2_2
                genos_DRBX_all.append(geno_DRBX_sorted)
                genos_DRB1_all.append(geno_DRB1)
                genos_DRBXDRB1_all.append(geno_phase2)

                # check for assoc rules
                # print ("Phase 2: DRB1 " + drb1_2_2dig + " DRBX " + drbx_loc1 + " DRB1 " + drb1_1_2dig +
                #         " DRBX " +  drbx_loc2)
                if ((DRBX_DRB1_assoc[drb1_2_2dig] == drbx_loc1) &
                    (DRBX_DRB1_assoc[drb1_1_2dig] == drbx_loc2)):
                    genos_DRBX_restricted.append(geno_DRBX_sorted)
                    genos_DRB1_restricted.append(geno_DRB1)
                    # print ("Accepted!")

        # check to see if any genotypes met association rules
        n_restricted_genos_DRBX = len(genos_DRBX_restricted)
        n_restricted_genos_DRB1 = len(genos_DRB1_restricted)
        # print (id n_restricted_genos_DRBX n_restricted_genos_DRB1)

        if (n_restricted_genos_DRB1 > 0):
            for geno in genos_DRBX_restricted:
                genos_id_DRBX[id][geno] = 1
            for geno in genos_DRB1_restricted:
                genos_id_DRB1[id][geno] = 1
        else:
            for geno in genos_DRBX_all:
                genos_id_DRBX[id][geno] = 1
            for geno in genos_DRB1_all:
                genos_id_DRB1[id][geno] = 1
            for geno in genos_DRBXDRB1_all:
                restricted_DRBX_DRB1_file.write("Restricted geno: " + id + " " + geno + "\n")

    restricted_DRBX_DRB1_file.close()

    # wipe out old glstring array
    glstring.clear()

    # generate new GLIDS for DRBX and DRB1 based on restricted list
    GLstrings_DRBX = {} # GLstrings made from block
    GLstrings_DRB1 = {} # GLstrings made from block
    GLID_DRBX = {} # GLID make from block for each donor
    GLID_DRB1 = {} # GLID make from block for each donor
    GLID_counter_DRBX = 3000000
    GLID_counter_DRB1 = 4000000
    glstring_DRBX_new = {}
    glstring_DRB1_new = {}
    # create new GLIDs for DRBX genotype lists
    for id in genos_id_DRBX:

        # print STDERR "$id\n";

        # sorts to limit unique GLIDs created
        GL = {}
        for geno in genos_id_DRBX[id]:
            GL[geno] = 1
        glstring_DRBX = "|".join(sorted(GL.keys()))

        # if GLstring is new, create new GLID
        if glstring_DRBX not in GLstrings_DRBX:    
            glstring_DRBX_new[GLID_counter_DRBX] = glstring_DRBX
            GLID_DRBX[id] = GLID_counter_DRBX
            GLstrings_DRBX[glstring_DRBX] = GLID_counter_DRBX
            GLID_counter_DRBX += 1
        else:
            GLID_DRBX[id] = GLstrings_DRBX[glstring_DRBX]

    # create new GLIDs for DRB1 genotype lists
    for id in genos_id_DRB1:

        # sort keys to limit unique GLIDs created
        GL = {}
        for geno in genos_id_DRB1[id]:
            GL[geno] = 1
        glstring_DRB1 = "|".join(sorted(GL.keys()))

        # if GLstring is new, create new GLID
        if glstring_DRB1 not in GLstrings_DRB1:    
            glstring_DRB1_new[GLID_counter_DRB1] = glstring_DRB1
            GLID_DRB1[id] = GLID_counter_DRB1
            GLstrings_DRB1[glstring_DRB1] = GLID_counter_DRB1
            GLID_counter_DRB1 += 1
        else:
            GLID_DRB1[id] = GLstrings_DRB1[glstring_DRB1]
        
    # replace pull ID dictionary with new GLIDs
    ID.clear()
    GLIDlist_count.clear()
    for id in GLID_DRB1:
        glid_drbx = GLID_DRBX[id]
        glid_drb1 = GLID_DRB1[id]
        glidlist = str(glid_drbx) + ":" + str(glid_drb1)
        ID[id] = glidlist
        glstring[glid_drb1] = glstring_DRB1_new[glid_drb1]
        glstring[glid_drbx] = glstring_DRBX_new[glid_drbx]
        GLIDlist_count[glidlist] += 1

# end if DRBX~DRB1

# DQA1~DQB1 two-locus block
if ((loc1 == "DQA1") & (loc2 == "DQB1")):

    DQA1_DQB1_assoc = {}
    DQA1_DQB1_assoc_filename = "./cfg/DQA1_DQB1_assoc.cfg"
    DQA1_DQB1_assoc_file = open (DQA1_DQB1_assoc_filename,'r')
    for line in DQA1_DQB1_assoc_file:
        line = line.rstrip('\r\n')
        (dqa1,dqb1) = line.split(',')
        if (dqb1 == "DQB1"):
            continue
        dq_haplo = dqa1 + "~" + dqb1
        DQA1_DQB1_assoc[dq_haplo] = 1
        # print STDERR "$dpa1 $dpb1 $DPA1_DPB1_assoc{$dpb1}\n";

    DQA1_DQB1_assoc_file.close()


    GLID_DQA1 = {}
    GLID_DQB1 = {}
    for id in ID:
        glidlist = ID[id]
        (glid_DQA1,glid_DQB1) = glidlist.split(":")
        if (glid_DQB1 != "999999999"):
            GLID_DQA1[id] = glid_DQA1
            GLID_DQB1[id] = glid_DQB1

    # replace pull ID dictionary with subset of GLIDs
    ID.clear()
    GLIDlist_count.clear()
    for id in GLID_DQB1:
        glid_dqa1 = GLID_DQA1[id]
        glid_dqb1 = GLID_DQB1[id]
        glidlist = str(glid_dqa1) + ":" + str(glid_dqb1)
        ID[id] = glidlist
        GLIDlist_count[glidlist] += 1

# DPA1~DPB1 two-locus block
if ((loc1 == "DPA1") & (loc2 == "DPB1")):

    DPA1_DPB1_assoc = {}
    DPA1_DPB1_assoc_filename = "./cfg/DPA1_DPB1_AA_assoc.cfg"
    DPA1_DPB1_assoc_file = open (DPA1_DPB1_assoc_filename,'r')
    for line in DPA1_DPB1_assoc_file:
        line = line.rstrip('\r\n')
        (dpa1,dpb1,dpa1_motif,dpb1_motif,stability) = line.split(',')
        if (dpb1 == "DPB1"):
            continue
        dp_haplo = dpa1 + "~" + dpb1
        DPA1_DPB1_assoc[dp_haplo] = stability
        # print STDERR "$dpa1 $dpb1 $DPA1_DPB1_assoc{$dpb1}\n";

    DPA1_DPB1_assoc_file.close()

    GLID_DPA1 = {}
    GLID_DPB1 = {}
    for id in ID:
        glidlist = ID[id]
        (glid_DPA1,glid_DPB1) = glidlist.split(":")
        if (glid_DPB1 != "777777777"):
            GLID_DPA1[id] = glid_DPA1
            GLID_DPB1[id] = glid_DPB1

    # replace pull ID dictionary with subset of GLIDs
    ID.clear()
    GLIDlist_count.clear()
    for id in GLID_DPB1:
        glid_dpa1 = GLID_DPA1[id]
        glid_dpb1 = GLID_DPB1[id]
        glidlist = str(glid_dpa1) + ":" + str(glid_dpb1)
        ID[id] = glidlist
        GLIDlist_count[glidlist] += 1

# print new GLID file
new_glid_filename = "./" + loc_folder_out +  "/glid." + dataset + "." + pop + ".txt"
new_glid_file = open(new_glid_filename,'w')
for glid in glstring:
    new_glid_file.write(str(glid) + ',' + glstring[glid] + "\n")
new_glid_file.close()

# print new donor pull file
new_pull_filename = "./" + loc_folder_out + "/pull." + dataset + "." + pop + ".txt"
new_pull_file = open(new_pull_filename,'w')
for id in ID:
    new_pull_file.write(id + "," + ID[id] + "\n")
new_pull_file.close()

# print EM input file
nemo_filename = "./" + loc_folder_out + "/nemo." + dataset + "." + pop + ".txt"
nemo_file = open(nemo_filename,'w')
for glidlist in GLIDlist_count:
    nemo_file.write (glidlist + "," + str(GLIDlist_count[glidlist]) + "\n")
nemo_file.close()

sys.stderr.write("--- %s seconds ---\n" % (time.time() - start_time))

exit (0)