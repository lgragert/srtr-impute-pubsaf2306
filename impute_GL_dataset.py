#!/usr/bin/env python
# impute_GL - impute genotype lists in NEMO pipeline

from collections import defaultdict
import sys
import copy
import time
import os
import gzip

print (sys.argv)

dataset = sys.argv[1]
pull_filename = sys.argv[2]
glid_filename = sys.argv[3]
HF_filename = sys.argv[4]
impute_filename = sys.argv[5]
impute_threshold = sys.argv[6]
max_genos = sys.argv[7]

sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)

start_time = time.time()

def convert_haplo_1F(haplo):
  alleles = haplo.split('~')
  alleles_1F = []
  for loctyp in alleles:
    AF_US[loctyp] += float(freq)
    if (loctyp == "DRBX*NNNN"):
      allele_1F = "DRBX*NN"
    else:
      (loc,typ) = loctyp.split('*')
      (typ_1F,_typ_2F) = typ.split(':')
      allele_1F = '*'.join([loc,typ_1F])

      # TODO - Find alternative to 1st field for DPB1 allele families, but use two-field for now
      if (loc == "DPB1"):
        allele_1F = loctyp
        
    alleles_1F.append(allele_1F)
  haplo_1F = '~'.join(alleles_1F)
  return (haplo_1F)




# load GLIDs

glstring = {} # stores GL String for each GLID
glid_file = open (glid_filename,'r')
for line in glid_file:
  line = line.rstrip('\r\n')
  (glid,gl) = line.split(',')
  glstring[glid] = gl

glid_file.close()

# get IDs from pull file
genos = {}
genos_ID = defaultdict(list)
GL = defaultdict(lambda: defaultdict(dict))
pull_file = open(pull_filename,'r')
for line in pull_file:
  line = line.rstrip('\r\n')
  (id, glidlist) = line.split(',')
  (glid1,glid2) = glidlist.split(':')
  genos[id] = glidlist
  genos_ID[glidlist].append(id) 
  gl1 = glstring[glid1].split('|') 
  for geno in gl1:
    if ("1" not in GL):
      GL["1"] = {}
    if (glid1 not in GL["1"]):
      GL["1"][glid1] = {}
    if (geno not in GL["1"][glid1]):
      GL["1"][glid1][geno] = 0
    GL["1"][glid1][geno] += 1

  gl2 = glstring[glid2].split('|')
  for geno in gl2:
    if ("2" not in GL):
      GL["2"] = {}
    if (glid2 not in GL["2"]):
      GL["2"][glid2] = {}
    if (geno not in GL["2"][glid2]):
      GL["2"][glid2][geno] = 0
    GL["2"][glid2][geno] += 1

pull_file.close()

# load HF
HF = {}
HF_1F = defaultdict(float) # first field haplotype freqs
AF = defaultdict(float)  # allele freqs for SHF
ndonors = 0
HF_file = open(HF_filename,'r')
for line in HF_file:
  (haplo,count,freq) = line.split(',')
  if (haplo == "Haplo"):
    continue
  if (float(freq) < 1E-10): # skip very rare haplotypes - might screw up imputation
    continue
  HF[haplo] = float(freq)
  ndonors += float(count) / 2

  # load first field haplotype freqs
  alleles = haplo.split('~')
  alleles_1F = []
  for loctyp in alleles:
    AF[loctyp] += float(freq)
    if (loctyp == "DRBX*NNNN"):
      allele_1F = "DRBX*NN"
    else:
      (loc,typ) = loctyp.split('*')
      (typ_1F,typ_2F) = loctyp.split(':')
      allele_1F = '*'.join([loc,typ_1F])
    alleles_1F.append(allele_1F)
  haplo_1F = '*'.join(alleles_1F)
  HF_1F[haplo_1F] += float(freq)
  
HF_file.close()

# parse out freqs directory and population name
(freqs_dir_list) = HF_filename.split('/')
freqs_dir = freqs_dir_list[1] # extract frequency directory
(freqs_file_list) = freqs_dir_list[2].split('.')
pop = freqs_file_list[1] # extract population name

# load US haplotypes for 2nd pass
HF_US = {}
HF_US_1F = defaultdict(float) # first field haplotype freqs
AF_US = defaultdict(float)  # allele freqs for SHF
HF_US_file = open("./" + freqs_dir + "/freqs.US.csv",'r')
for line in HF_US_file:
  (haplo,count,freq) = line.split(',')
  if (haplo == "Haplo"):
    continue
  HF_US[haplo] = float(freq)
  if (float(freq) < 1E-10):  # skip very rare haplotypes - might screw up imputation
    continue

  # load first field haplotype freqs
  alleles = haplo.split('~')
  alleles_1F = []
  for loctyp in alleles:
    AF_US[loctyp] += float(freq)
    if (loctyp == "DRBX*NNNN"):
      allele_1F = "DRBX*NN"
    else:
      (loc,typ) = loctyp.split('*')
      (typ_1F,typ_2F) = typ.split(':')
      allele_1F = '*'.join([loc,typ_1F])

      # TODO - Find alternative to 1st field for DPB1 allele families, but use two-field for now
      if (loc == "DPB1"):
        allele_1F = loctyp
        
    alleles_1F.append(allele_1F)
  haplo_1F = '~'.join(alleles_1F)
  # print (haplo_1F)
  HF_US_1F[haplo_1F] += float(freq)
  
HF_US_file.close()


# TEMP load CAU haplotypes for 3rd pass
HF_CAU = {}
HF_CAU_1F = defaultdict(float) # first field haplotype freqs
AF_US = defaultdict(float)  # allele freqs for SHF
HF_CAU_file = open("./" + freqs_dir + "/freqs.CAU.csv",'r')
for line in HF_CAU_file:
  (haplo,count,freq) = line.split(',')
  if (haplo == "Haplo"):
    continue
  if (float(freq) < 1E-10):  # skip very rare haplotypes - might screw up imputation
    continue
  HF_CAU[haplo] = float(freq)

  # load first field haplotype freqs
  alleles = haplo.split('~')
  alleles_1F = []
  for loctyp in alleles:
    AF_US[loctyp] += float(freq)
    if (loctyp == "DRBX*NNNN"):
      allele_1F = "DRBX*NN"
    else:
      (loc,typ) = loctyp.split('*')
      (typ_1F,typ_2F) = typ.split(':')
      allele_1F = '*'.join([loc,typ_1F])

      # TODO - Find alternative to 1st field for DPB1 allele families, but use two-field for now
      if (loc == "DPB1"):
        allele_1F = loctyp
        
    alleles_1F.append(allele_1F)
  haplo_1F = '~'.join(alleles_1F)
  # print (haplo_1F)
  HF_CAU_1F[haplo_1F] += float(freq)
  
HF_CAU_file.close()


# make HF distribution incremented by HF of cases that did not impute
HF_cases = HF

sys.stderr.write("Imputation Input File: " + pull_filename + "\n")
sys.stderr.write("Load Time: %s seconds\n" % (time.time() - start_time))

# imputation output
glidlist_increment = {} # glidlists that don't impute on first pass get incremented to freqs
impute_file = gzip.open (impute_filename,'wt')
for glidlist in genos_ID:

  # geno = genos[id]
  (glid1,glid2) = glidlist.split(":")
  haplo_pair_freq = defaultdict(float)
  haplo_pair_total = 0
  for geno1 in GL["1"][glid1]:

    (a1,a2) = geno1.split("+")

    for geno2 in GL["2"][glid2]:

      (b1,b2) = geno2.split("+")

      hap_11 = a1 + "~" + b1
      hap_12 = a1 + "~" + b2
      hap_21 = a2 + "~" + b1
      hap_22 = a2 + "~" + b2

      # DRB1DQB1 and DQA1:
      if (b1[0:4] == "DQA1") & (a1[0:4] == "DRB1"):
        dqa1_1 = b1
        dqa1_2 = b2
        (drb1_1,dqb1_1) = a1.split("~")
        (drb1_2,dqb1_2) = a2.split("~") 
        hap_11 = '~'.join([drb1_1,dqa1_1,dqb1_1])
        hap_12 = '~'.join([drb1_1,dqa1_2,dqb1_1])
        hap_21 = '~'.join([drb1_2,dqa1_1,dqb1_2])
        hap_22 = '~'.join([drb1_2,dqa1_2,dqb1_2])

      # DRBXDRB1DQB1 and DQA1
      if ((b1[0:4] == "DQA1") & (a1[0:4] in ["DRB3","DRB4","DRB5","DRBX"])):
        dqa1_1 = b1
        dqa1_2 = b2
        (drbx_1,drb1_1,dqb1_1) = a1.split("~")
        (drbx_2,drb1_2,dqb1_2) = a2.split("~") 
        hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1])
        hap_12 = '~'.join([drbx_1,drb1_1,dqa1_2,dqb1_1])
        hap_21 = '~'.join([drbx_2,drb1_2,dqa1_1,dqb1_2])
        hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2])

      # DRBXDRB1DQA1DQB1DPB1 and DPA1
      if (b1[0:4] == "DPA1"):
        dpa1_1 = b1
        dpa1_2 = b2
        (drbx_1,drb1_1,dqa1_1,dqb1_1,dpb1_1) = a1.split("~")
        (drbx_2,drb1_2,dqa1_2,dqb1_2,dpb1_2) = a2.split("~") 
        hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_1,dpb1_1])
        hap_12 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_2,dpb1_1])
        hap_21 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_1,dpb1_2])
        hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_2,dpb1_2])

      # hap_11 = "~".join([a1,b1])
      # hap_12 = "~".join([a1,b2])
      # hap_21 = "~".join([a2,b1])
      # hap_22 = "~".join([a2,b2])

      # FIRST PASS - Set haplotype freq to zero for missing haplotypes and see if they impute
      if (hap_11 not in HF):
        HF[hap_11] = 0
      if (hap_12 not in HF):
        HF[hap_12] = 0
      if (hap_21 not in HF):
        HF[hap_21] = 0
      if (hap_22 not in HF):
        HF[hap_22] = 0  

      if (hap_11 > hap_22):
        (hap_11, hap_22) = (hap_22, hap_11)

      if (hap_12 > hap_21):
        (hap_12, hap_21) = (hap_21, hap_12)        

      # two possible pairs of genotypes
      if (hap_11 != hap_22):
        freq = 2 * HF[hap_11] * HF[hap_22]
      else:
        freq = HF[hap_11] ** 2

      if (freq != 0):
        haplo_pair_freq[hap_11 + "+" + hap_22] += freq
        haplo_pair_total += freq

      if (hap_12 != hap_21):
        freq = 2 * HF[hap_12] * HF[hap_21]
      else:
        freq = HF[hap_12] ** 2

      if (freq != 0):
        haplo_pair_freq[hap_12 + "+" + hap_21] += freq
      # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))
        haplo_pair_total += freq


  # SECOND PASS - Run imputation using CAU freqs
  if (len(haplo_pair_freq) == 0):
    sys.stderr.write("Second Pass CAU Freqs: " + genos_ID[glidlist][0] + ": " + glidlist + "  " + glstring[glid1] + " " + glstring[glid2] + "\n")
    for geno1 in GL["1"][glid1]:

      (a1,a2) = geno1.split("+")

      for geno2 in GL["2"][glid2]:

        (b1,b2) = geno2.split("+")

        hap_11 = a1 + "~" + b1
        hap_12 = a1 + "~" + b2
        hap_21 = a2 + "~" + b1
        hap_22 = a2 + "~" + b2

        # DRB1DQB1 and DQA1:
        if (b1[0:4] == "DQA1") & (a1[0:4] == "DRB1"):
          dqa1_1 = b1
          dqa1_2 = b2
          (drb1_1,dqb1_1) = a1.split("~")
          (drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQB1 and DQA1
        if ((b1[0:4] == "DQA1") & (a1[0:4] in ["DRB3","DRB4","DRB5","DRBX"])):
          dqa1_1 = b1
          dqa1_2 = b2
          (drbx_1,drb1_1,dqb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQA1DQB1DPB1 and DPA1
        if (b1[0:4] == "DPA1"):
          dpa1_1 = b1
          dpa1_2 = b2
          (drbx_1,drb1_1,dqa1_1,dqb1_1,dpb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqa1_2,dqb1_2,dpb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_1,dpb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_2,dpb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_1,dpb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_2,dpb1_2])

        # SECOND PASS - Set US haplotype freq to zero for missing haplotypes and see if they impute
        if (hap_11 not in HF_CAU):
          HF_CAU[hap_11] = 0
        if (hap_12 not in HF_CAU):
          HF_CAU[hap_12] = 0
        if (hap_21 not in HF_CAU):
          HF_CAU[hap_21] = 0
        if (hap_22 not in HF_CAU):
          HF_CAU[hap_22] = 0  

        # hap_11 = "~".join([a1,b1])
        # hap_12 = "~".join([a1,b2])
        # hap_21 = "~".join([a2,b1])
        # hap_22 = "~".join([a2,b2])

        if (hap_11 > hap_22):
          (hap_11, hap_22) = (hap_22, hap_11)

        if (hap_12 > hap_21):
          (hap_12, hap_21) = (hap_21, hap_12)        

        # two possible pairs of genotypes
        if (hap_11 != hap_22):
          freq = 2 * HF_CAU[hap_11] * HF_CAU[hap_22]
        else:
          freq = HF_CAU[hap_11] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_11 + "+" + hap_22] += freq
          haplo_pair_total += freq

        if (hap_12 != hap_21):
          freq = 2 * HF_CAU[hap_12] * HF_CAU[hap_21]
        else:
          freq = HF_CAU[hap_12] ** 2

        # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))

        if (freq != 0):
          haplo_pair_freq[hap_12 + "+" + hap_21] += freq
        # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))
          haplo_pair_total += freq

  # SECOND PASS - Run imputation using US freqs
  if (len(haplo_pair_freq) == 0):
    sys.stderr.write("Second Pass US Freqs: " + genos_ID[glidlist][0] + ": " + glidlist + "  " + glstring[glid1] + " " + glstring[glid2] + "\n")
    for geno1 in GL["1"][glid1]:

      (a1,a2) = geno1.split("+")

      for geno2 in GL["2"][glid2]:

        (b1,b2) = geno2.split("+")

        hap_11 = a1 + "~" + b1
        hap_12 = a1 + "~" + b2
        hap_21 = a2 + "~" + b1
        hap_22 = a2 + "~" + b2

        # DRB1DQB1 and DQA1:
        if (b1[0:4] == "DQA1") & (a1[0:4] == "DRB1"):
          dqa1_1 = b1
          dqa1_2 = b2
          (drb1_1,dqb1_1) = a1.split("~")
          (drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQB1 and DQA1
        if ((b1[0:4] == "DQA1") & (a1[0:4] in ["DRB3","DRB4","DRB5","DRBX"])):
          dqa1_1 = b1
          dqa1_2 = b2
          (drbx_1,drb1_1,dqb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQA1DQB1DPB1 and DPA1
        if (b1[0:4] == "DPA1"):
          dpa1_1 = b1
          dpa1_2 = b2
          (drbx_1,drb1_1,dqa1_1,dqb1_1,dpb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqa1_2,dqb1_2,dpb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_1,dpb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_2,dpb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_1,dpb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_2,dpb1_2])

        # SECOND PASS - Set US haplotype freq to zero for missing haplotypes and see if they impute
        if (hap_11 not in HF_CAU):
          HF_CAU[hap_11] = 0
        if (hap_12 not in HF_CAU):
          HF_CAU[hap_12] = 0
        if (hap_21 not in HF_CAU):
          HF_CAU[hap_21] = 0
        if (hap_22 not in HF_CAU):
          HF_CAU[hap_22] = 0  

        # hap_11 = "~".join([a1,b1])
        # hap_12 = "~".join([a1,b2])
        # hap_21 = "~".join([a2,b1])
        # hap_22 = "~".join([a2,b2])

        if (hap_11 > hap_22):
          (hap_11, hap_22) = (hap_22, hap_11)

        if (hap_12 > hap_21):
          (hap_12, hap_21) = (hap_21, hap_12)        

        # two possible pairs of genotypes
        if (hap_11 != hap_22):
          freq = 2 * HF_CAU[hap_11] * HF_CAU[hap_22]
        else:
          freq = HF_CAU[hap_11] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_11 + "+" + hap_22] += freq
          haplo_pair_total += freq

        if (hap_12 != hap_21):
          freq = 2 * HF_CAU[hap_12] * HF_CAU[hap_21]
        else:
          freq = HF_CAU[hap_12] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_12 + "+" + hap_21] += freq
        # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))
          haplo_pair_total += freq


  # THIRD PASS - Run SHF imputation using first-field haplotype freqs
  if (len(haplo_pair_freq) == 0):
    sys.stderr.write("Third Pass - Two-Field HF: " + genos_ID[glidlist][0] + ": " + glidlist + "  " + glstring[glid1] + " " + glstring[glid2] + "\n")
    for geno1 in GL["1"][glid1]:

      (a1,a2) = geno1.split("+")

      for geno2 in GL["2"][glid2]:

        (b1,b2) = geno2.split("+")

        hap_11 = a1 + "~" + b1
        hap_12 = a1 + "~" + b2
        hap_21 = a2 + "~" + b1
        hap_22 = a2 + "~" + b2

        # DRB1DQB1 and DQA1:
        if (b1[0:4] == "DQA1") & (a1[0:4] == "DRB1"):
          dqa1_1 = b1
          dqa1_2 = b2
          (drb1_1,dqb1_1) = a1.split("~")
          (drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQB1 and DQA1
        if ((b1[0:4] == "DQA1") & (a1[0:4] in ["DRB3","DRB4","DRB5","DRBX"])):
          dqa1_1 = b1
          dqa1_2 = b2
          (drbx_1,drb1_1,dqb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQA1DQB1DPB1 and DPA1
        if (b1[0:4] == "DPA1"):
          dpa1_1 = b1
          dpa1_2 = b2
          (drbx_1,drb1_1,dqa1_1,dqb1_1,dpb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqa1_2,dqb1_2,dpb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_1,dpb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_2,dpb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_1,dpb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_2,dpb1_2])

        # THIRD PASS - Use first-field haplotypes as fallback - Doesn't distinguish among alleles
        if (hap_11 not in HF):
          hap_11_1F = convert_haplo_1F(hap_11)
          HF[hap_11] = HF_US_1F[hap_11_1F]
        if (hap_12 not in HF):
          hap_12_1F = convert_haplo_1F(hap_12)
          HF[hap_12] = HF_US_1F[hap_12_1F]
        if (hap_21 not in HF):
          hap_21_1F = convert_haplo_1F(hap_21)
          HF[hap_21] = HF_US_1F[hap_21_1F]
        if (hap_22 not in HF):
          hap_22_1F = convert_haplo_1F(hap_22)
          HF[hap_22] = HF_US_1F[hap_22_1F]

        if (hap_11 > hap_22):
          (hap_11, hap_22) = (hap_22, hap_11)

        if (hap_12 > hap_21):
          (hap_12, hap_21) = (hap_21, hap_12)        

        # two possible pairs of genotypes
        if (hap_11 != hap_22):
          freq = 2 * HF[hap_11] * HF[hap_22]
        else:
          freq = HF[hap_11] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_11 + "+" + hap_22] += freq
          haplo_pair_total += freq

        if (hap_12 != hap_21):
          freq = 2 * HF[hap_12] * HF[hap_21]
        else:
          freq = HF[hap_12] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_12 + "+" + hap_21] += freq
        # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))
          haplo_pair_total += freq


  # FOURTH PASS - Run SHF imputation using allele freqs for cases where no first-field data
  if (len(haplo_pair_freq) == 0):
    sys.stderr.write("Fourth Pass - Allele Freq: " + genos_ID[glidlist][0] + ": " + glidlist + "  " + glstring[glid1] + " " + glstring[glid2] + "\n")
    for geno1 in GL["1"][glid1]:

      (a1,a2) = geno1.split("+")

      for geno2 in GL["2"][glid2]:

        (b1,b2) = geno2.split("+")

        hap_11 = a1 + "~" + b1
        hap_12 = a1 + "~" + b2
        hap_21 = a2 + "~" + b1
        hap_22 = a2 + "~" + b2

        # DRB1DQB1 and DQA1:
        if (b1[0:4] == "DQA1") & (a1[0:4] == "DRB1"):
          dqa1_1 = b1
          dqa1_2 = b2
          (drb1_1,dqb1_1) = a1.split("~")
          (drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQB1 and DQA1
        if ((b1[0:4] == "DQA1") & (a1[0:4] in ["DRB3","DRB4","DRB5","DRBX"])):
          dqa1_1 = b1
          dqa1_2 = b2
          (drbx_1,drb1_1,dqb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQA1DQB1DPB1 and DPA1
        if (b1[0:4] == "DPA1"):
          dpa1_1 = b1
          dpa1_2 = b2
          (drbx_1,drb1_1,dqa1_1,dqb1_1,dpb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqa1_2,dqb1_2,dpb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_1,dpb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_2,dpb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_1,dpb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_2,dpb1_2])

        if (a1 not in AF):
          AF[a1] = 0
        if (a2 not in AF):
          AF[a2] = 0
        if (b1 not in AF):
          AF[b1] = 0
        if (b2 not in AF):
          AF[b2] = 0

        # FOURTH PASS - SHF for haplotypes - Use allele frequencies as a fallback
        if (hap_11 not in HF):
          HF[hap_11] = AF[a1] * AF[b1]
        if (hap_12 not in HF):
          HF[hap_12] = AF[a1] * AF[b2]
        if (hap_21 not in HF):
          HF[hap_21] = AF[a2] * AF[b1]
        if (hap_22 not in HF):
          HF[hap_22] = AF[a2] * AF[b2]

        if (hap_11 > hap_22):
          (hap_11, hap_22) = (hap_22, hap_11)

        if (hap_12 > hap_21):
          (hap_12, hap_21) = (hap_21, hap_12)        

        # two possible pairs of genotypes
        if (hap_11 != hap_22):
          freq = 2 * HF[hap_11] * HF[hap_22]
        else:
          freq = HF[hap_11] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_11 + "+" + hap_22] += freq
          haplo_pair_total += freq

        if (hap_12 != hap_21):
          freq = 2 * HF[hap_12] * HF[hap_21]
        else:
          freq = HF[hap_12] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_12 + "+" + hap_21] += freq
        # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))
          haplo_pair_total += freq

  # FIFTH AND FINAL PASS - Initialize allele freq at low probability when is missing completely
  if (len(haplo_pair_freq) == 0):
    sys.stderr.write("Fifth Pass - New Alleles: " + genos_ID[glidlist][0] + ": " + glidlist + "  " + glstring[glid1] + " " + glstring[glid2] + "\n")
    for geno1 in GL["1"][glid1]:

      (a1,a2) = geno1.split("+")

      for geno2 in GL["2"][glid2]:

        (b1,b2) = geno2.split("+")

        hap_11 = a1 + "~" + b1
        hap_12 = a1 + "~" + b2
        hap_21 = a2 + "~" + b1
        hap_22 = a2 + "~" + b2

        # DRB1DQB1 and DQA1:
        if (b1[0:4] == "DQA1") & (a1[0:4] == "DRB1"):
          dqa1_1 = b1
          dqa1_2 = b2
          (drb1_1,dqb1_1) = a1.split("~")
          (drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQB1 and DQA1
        if ((b1[0:4] == "DQA1") & (a1[0:4] in ["DRB3","DRB4","DRB5","DRBX"])):
          dqa1_1 = b1
          dqa1_2 = b2
          (drbx_1,drb1_1,dqb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_2,dqb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_1,dqb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2])

        # DRBXDRB1DQA1DQB1DPB1 and DPA1
        if (b1[0:4] == "DPA1"):
          dpa1_1 = b1
          dpa1_2 = b2
          (drbx_1,drb1_1,dqa1_1,dqb1_1,dpb1_1) = a1.split("~")
          (drbx_2,drb1_2,dqa1_2,dqb1_2,dpb1_2) = a2.split("~") 
          hap_11 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_1,dpb1_1])
          hap_12 = '~'.join([drbx_1,drb1_1,dqa1_1,dqb1_1,dpa1_2,dpb1_1])
          hap_21 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_1,dpb1_2])
          hap_22 = '~'.join([drbx_2,drb1_2,dqa1_2,dqb1_2,dpa1_2,dpb1_2])

        # FOURTH PASS - default allele frequency
        if (AF[a1] == 0):
          AF[a1] = 0.00000001
        if (AF[a2] == 0):
          AF[a2] = 0.00000001
        if (AF[b1] == 0):
          AF[b1] = 0.00000001
        if (AF[b2] == 0):
          AF[b2] = 0.00000001
        # TODO - Use two-field freqs as a multiplier except for DPB1
        if (HF[hap_11] == 0):
          HF[hap_11] = AF[a1] * AF[b1]
        if (HF[hap_12] == 0):
          HF[hap_12] = AF[a1] * AF[b2]
        if (HF[hap_21] == 0):
          HF[hap_21] = AF[a2] * AF[b1]
        if (HF[hap_22] == 0):
          HF[hap_22] = AF[a2] * AF[b2]

        # hap_11 = "~".join([a1,b1])
        # hap_12 = "~".join([a1,b2])
        # hap_21 = "~".join([a2,b1])
        # hap_22 = "~".join([a2,b2])

        if (hap_11 > hap_22):
          (hap_11, hap_22) = (hap_22, hap_11)

        if (hap_12 > hap_21):
          (hap_12, hap_21) = (hap_21, hap_12)        

        # two possible pairs of genotypes
        if (hap_11 != hap_22):
          freq = 2 * HF[hap_11] * HF[hap_22]
        else:
          freq = HF[hap_11] ** 2

        # print (id + " " + geno1 + " " + geno2 + " " + hap_11 + " " + hap_12 + " " + hap_21 + " " + hap_22 + " " + str(freq))

        if (freq != 0):
          haplo_pair_freq[hap_11 + "+" + hap_22] += freq
          haplo_pair_total += freq

        if (hap_12 != hap_21):
          freq = 2 * HF[hap_12] * HF[hap_21]
        else:
          freq = HF[hap_12] ** 2

        if (freq != 0):
          haplo_pair_freq[hap_12 + "+" + hap_21] += freq
          haplo_pair_total += freq

  # Final check for failed imputation
  if (len(haplo_pair_freq) == 0):
    print ("Did not impute: " + id)

  geno_count = 0
  prob_cumul = 0


  impute_out_array = []
  for hap_pair, GF in sorted(haplo_pair_freq.items(), key=lambda item: item[1], reverse=True):
    if (prob_cumul > float(impute_threshold)):
      # print ("Impute threshold reached: " + str(prob_cumul))
      continue
    if (geno_count > int(max_genos)):
      continue
    prob = GF / haplo_pair_total
    # print (hap_pair)
    (hap1, hap2) = hap_pair.split('+')
    prob_cumul += prob
    geno_count += 1
    # print (id + ',' + geno + "," + str(geno_count) + "," + hap_pair + "," + str(GF) + "," + hap1 + "," 
    #       + str(HF[hap1]) + "," + hap2 + "," + str(HF[hap2]) + "," + str(prob) + "," + str(prob_cumul))
    impute_out_array.append(str(geno_count) + "," + hap1 + "," + hap2 + "," + str(prob))

  for id in genos_ID[glidlist]:
    for impute_out_line in impute_out_array:
      impute_file.write(id + "," + impute_out_line + "\n")

  # HF
  if glidlist in glidlist_increment:
    for impute_out_line in impute_out_array:
      (geno_count_impute,hap1_impute,hap2_impute,prob_impute) = impute_out_line.split(',')
      HF_cases[hap1_impute] += geno_count_impute / ndonors * prob_impute
      HF_cases[hap2_impute] += geno_count_impute / ndonors * prob_impute

impute_file.close()

sys.stderr.write("--- %s seconds ---\n" % (time.time() - start_time))

exit (0)

# write out new haplotype freqs with incremented freqs
HF_dataset_filename = "./" + freqs_dir + "/freqs." + dataset + "." + pop + ".csv"
HF_dataset_file = open(HF_cases_filename,'w')
HF_dataset_file.write("Haplo,Count,Freq\n")

HF_cases_lines = []
for hap in HF_cases:
  freq = HF_cases[hap]
  count = HF_cases[hap]/ (2 * ndonors)
  HF_cases_lines.append(hap + ',' + str(count) + ',' + str(freq) + "\n")

HF_cases_file.writelines(HF_cases_lines)

HF_cases_file.close()


