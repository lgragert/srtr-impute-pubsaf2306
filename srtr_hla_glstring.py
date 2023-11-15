#!/usr/bin/env python
import pandas
import numpy
import re
import requests
import pyard

# initialize pyARD and roll up typing to 'g' group for imputation
# XX codes are always in latest IMGT/HLA allele list
ard = pyard.ARD()

# global XX code GL hash
XX_hash = {}
dna2ser = {}
seroGL = {}

def expand_XX(ac, XX_hash):
	# uses hashing to limit calls to MAC service 
	if (ac in XX_hash):
		gl = XX_hash[ac]
		# print ("Hash Lookup")
	else:
		url = "https://hml.nmdp.org/mac/api/decode?typing="
		response = requests.get(url + ac)
		gl = response.text
		XX_hash[ac] = gl
		# print ("Service Call")

	gl_ars = ard.redux_gl(gl,'lgx')
	return(gl_ars)

def expand_AC(ac):
	# print ("AC " + ac)
	gl_ars = ard.redux_gl(ac,'lg')
	return(gl_ars)

def loc_gl (loc, typ1, typ2, dna1, dna2, XX_hash):
	if (typ2 == ""):
		typ2 = typ1
		dna2 = dna1
	if (typ1 == ""):
		typ1 = typ2
		dna1 = dna2
	if ((typ1 == "") | (typ2 == "")):
		# if (loc in ["A","B","DRB1"]):
		# 	print ("Either no typ1 or typ2 in loc_gl: " + str(getattr(row, 'PX_ID')) + " Typ1: " + typ1 + " Typ2: " + typ2 + "^")
		return "SKIP"

	# print (loc + " TYP: " + typ1 + " " + typ2 + " DNA: " + dna)	
	
	al1 = ""
	al2 = ""

	# C is DNA only
	if (loc == "C"):
		dna1 = "Y"
		dna2 = "Y"
		if ((typ1 == "13:XX") | (typ2 == "13:XX")):
			# Cw*1301	Sequence shown to be in error (July 2002) - https://hla.alleles.org/alleles/deleted.html
			# print ("Deleted allele C*13:XX in HLA typing: " + str(getattr(row, 'PX_ID')))
			return "SKIP"

	# B*15:22 renamed to B*35:43
	if (loc == "B"):
		if (typ1 == "15:22"):
			typ1 = "35:43"
		if (typ2 == "15:22"):
			typ2 = "35:43"			


	# serology
	if (dna1 == "N"):
		loc_ser = loc
		if (loc == "DRB1"):
			loc_ser = "DR"
		if (loc == "DQB1"):
			loc_ser = "DQ"
		if (loc == "DPB1"):
			loc_ser = "DP"
			# print ("Legacy DP antigen-level typing  - not using for imputation:" + str(getattr(row, 'PX_ID')))
			return "SKIP"  # DP serology is inactive
		loctyp1 = loc_ser + typ1
		# print ("SERO " + loctyp1 + "" + loctyp2)
		al1 = seroGL[loctyp1]

	# dna typing
	elif (dna1 == "Y"):
		typ1_fields = typ1.split(':')
		subtype1 = typ1_fields[1]
		loctyp1 = '*'.join([loc,typ1])
		if (subtype1 == "XX"):
			al1 = expand_XX(loctyp1,XX_hash)
		else:
			al1 = loctyp1
			# ARD rollup of alleles in typing - Handle cases like DPB1*350:01 which aren't in freqs
			al1 = ard.redux_gl(al1,'lgx')

	if (dna2 == "N"):
		loc_ser = loc
		if (loc == "DRB1"):
			loc_ser = "DR"
		if (loc == "DQB1"):
			loc_ser = "DQ"
		if (loc == "DPB1"):
			loc_ser = "DP"
			# print ("Legacy DP antigen-level typing - not using for imputation: " + str(getattr(row, 'PX_ID')))
			return "SKIP"  # DP serology is inactive
		loctyp2 = loc_ser + typ2
		al2 = seroGL[loctyp2]

	# dna typing
	elif (dna2 == "Y"):
		typ2_fields = typ2.split(':')
		subtype2 = typ2_fields[1]
		loctyp2 = '*'.join([loc,typ2])
		if (subtype2 == "XX"):
			al2 = expand_XX(loctyp2,XX_hash)
		else:
			al2 = loctyp2
			# ARD rollup of alleles in typing - Handle cases like DPB1*350:01 which aren't in freqs
			al2 = ard.redux_gl(al2,'lgx')

	else:
		print ("Missing DNA or Serology Indicator for Locus Typing: " + str(getattr(row, 'PX_ID')))
		return "SKIP"

	gl = '+'.join([al1, al2])
	return gl


##############################################################################
# Function: loadSeroMap - load UNOS antigens to alleles
##############################################################################
def loadSeroWMDAMap (seroGL):

	unos_sero_filename = "OPTN_antigens_to_alleles_CPRA.txt"
	unos_sero_file = open(unos_sero_filename, 'r')

	seroGL["A203"] = "A*02:03"
	seroGL["A210"] = "A*02:10"
	seroGL["B703"] = "B*07:03"
	seroGL["B83"] = "B*83:01"

	for line in unos_sero_file:
		line = line.strip("\n\r")
		(antigen,alleles) = line.split("\t")
		seroGL[antigen] = alleles
		# print (seroGL[antigen])

	unos_sero_file.close()

	return (0)




# load file generated from HLA merge
tx_ki_hla_filename = "tx_ki_hla.csv"
tx_ki_hla = pandas.read_csv(tx_ki_hla_filename, index_col=None, encoding='Latin-1', low_memory=False)

# output SRTR pull/glid files
pull_filename = "pull.srtr.txt"
PULL = open(pull_filename, "w")

glid_filename = "glid.srtr.txt"
GLID = open(glid_filename, "w")


loadSeroWMDAMap(seroGL)

# PULL.write ("SUBJECT_IID,ROLLUP_RACE_CDE,glid_a,glid_b,glid_c,glid_drb1,glid_drbx,glid_dqa1,glid_dqb1,glid_dpa1,glid_dpb1\n")
# GLID.write ("GLID,GLSTRING\n")

glid_dict = {}
glid_index = 0

for row in tx_ki_hla.itertuples():

	# handle unknown and missing race
	DON_RACE = str(getattr(row,'DON_RACE'))
	CAN_RACE = str(getattr(row,'CAN_RACE'))
	if (CAN_RACE == "1024: Unknown (for Donor Referral only)"):
		CAN_RACE = "UNK"
	if (DON_RACE == "1024: Unknown (for Donor Referral only)"):
		DON_RACE = "UNK"	
	if (CAN_RACE == "Missing"):
		CAN_RACE = "UNK"
	if (DON_RACE == "Missing"):
		DON_RACE = "UNK"	

	# handle nan
	DON_A1 = str(getattr(row,'DON_A1'))
	DON_A2 = str(getattr(row,'DON_A2'))
	DON_C1 = str(getattr(row,'DON_C1'))
	DON_C2 = str(getattr(row,'DON_C2'))
	DON_B1 = str(getattr(row,'DON_B1'))
	DON_B2 = str(getattr(row,'DON_B2'))
	DON_DR1 = str(getattr(row,'DON_DR1'))
	DON_DR2 = str(getattr(row,'DON_DR2'))
	DON_DQ1 = str(getattr(row,'DON_DQ1'))
	DON_DQ2 = str(getattr(row,'DON_DQ2'))
	DON_DP1 = str(getattr(row,'DON_DP1'))
	DON_DP2 = str(getattr(row,'DON_DP2'))
	DON_DR51 = str(getattr(row,'DON_DR51')) # Negative or Positive until recently
	DON_DR52 = str(getattr(row,'DON_DR52')) # Negative or Positive until recently
	DON_DR53 = str(getattr(row,'DON_DR53')) # Negative or Positive until recently
	if (DON_A1 == "nan"):
		DON_A1 = ""
	if (DON_A2 == "nan"):
		DON_A2 = ""
	if (DON_C1 == "nan"):
		DON_C1 = ""
	if (DON_C2 == "nan"):
		DON_C2 = ""
	if (DON_B1 == "nan"):
		DON_B1 = ""
	if (DON_B2 == "nan"):
		DON_B2 = ""
	if (DON_DR1 == "nan"):
		DON_DR1 = ""
	if (DON_DR2 == "nan"):
		DON_DR2 = ""
	if (DON_DQ1 == "nan"):
		DON_DQ1 = ""
	if (DON_DQ2 == "nan"):
		DON_DQ2 = ""
	if (DON_DP1 == "nan"):
		DON_DP1 = ""
	if (DON_DP2 == "nan"):
		DON_DP2 = ""
	if (DON_DR51 == "nan"):
		DON_DR51 = ""
	if (DON_DR52 == "nan"):
		DON_DR52 = ""	
	if (DON_DR53 == "nan"):
		DON_DR53 = ""

	DON_A1_DNA_TYPING_IND = "N"
	DON_A2_DNA_TYPING_IND = "N"
	DON_C1_DNA_TYPING_IND = "N"
	DON_C2_DNA_TYPING_IND = "N"
	DON_B1_DNA_TYPING_IND = "N"
	DON_B2_DNA_TYPING_IND = "N"
	DON_DR1_DNA_TYPING_IND = "N"
	DON_DR2_DNA_TYPING_IND = "N"
	DON_DQ1_DNA_TYPING_IND = "N"
	DON_DQ2_DNA_TYPING_IND = "N"
	DON_DP1_DNA_TYPING_IND = "N"
	DON_DP2_DNA_TYPING_IND = "N"

	if (re.search(":",DON_A1)):
		DON_A1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_A2)):
		DON_A2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_C1)):
		DON_C1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_C2)):
		DON_C2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_B1)):
		DON_B1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_B2)):
		DON_B2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DR1)):
		DON_DR1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DR2)):
		DON_DR2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DQ1)):
		DON_DQ1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DQ2)):
		DON_DQ2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DP1)):
		DON_DP1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DP2)):
		DON_DP2_DNA_TYPING_IND = "Y"

	# handle nan
	REC_A1 = str(getattr(row,'REC_A1'))
	REC_A2 = str(getattr(row,'REC_A2'))
	REC_C1 = str(getattr(row,'REC_CW1'))
	REC_C2 = str(getattr(row,'REC_CW2'))
	REC_B1 = str(getattr(row,'REC_B1'))
	REC_B2 = str(getattr(row,'REC_B2'))
	REC_DR1 = str(getattr(row,'REC_DR1'))
	REC_DR2 = str(getattr(row,'REC_DR2'))
	REC_DQ1 = str(getattr(row,'REC_DQW1'))
	REC_DQ2 = str(getattr(row,'REC_DQW2'))
	REC_DP1 = str(getattr(row,'REC_DPW1'))
	REC_DP2 = str(getattr(row,'REC_DPW2'))
	REC_DR51 = str(getattr(row,'REC_DRW51')) # Negative or Positive until recently
	REC_DR52 = str(getattr(row,'REC_DRW52')) # Negative or Positive until recently
	REC_DR53 = str(getattr(row,'REC_DRW53')) # Negative or Positive until recently


	if (REC_A1 == "nan"):
		REC_A1 = ""
	if (REC_A2 == "nan"):
		REC_A2 = ""
	if (REC_C1 == "nan"):
		REC_C1 = ""
	if (REC_C2 == "nan"):
		REC_C2 = ""
	if (REC_B1 == "nan"):
		REC_B1 = ""
	if (REC_B2 == "nan"):
		REC_B2 = ""
	if (REC_DR1 == "nan"):
		REC_DR1 = ""
	if (REC_DR2 == "nan"):
		REC_DR2 = ""
	if (REC_DQ1 == "nan"):
		REC_DQ1 = ""
	if (REC_DQ2 == "nan"):
		REC_DQ2 = ""
	if (REC_DP1 == "nan"):
		REC_DP1 = ""
	if (REC_DP2 == "nan"):
		REC_DP2 = ""
	if (REC_DR51 == "nan"):
		REC_DR51 = ""
	if (REC_DR52 == "nan"):
		REC_DR52 = ""	
	if (REC_DR53 == "nan"):
		REC_DR53 = ""

	REC_A1_DNA_TYPING_IND = "N"
	REC_A2_DNA_TYPING_IND = "N"
	REC_C1_DNA_TYPING_IND = "N"
	REC_C2_DNA_TYPING_IND = "N"
	REC_B1_DNA_TYPING_IND = "N"
	REC_B2_DNA_TYPING_IND = "N"
	REC_DR1_DNA_TYPING_IND = "N"
	REC_DR2_DNA_TYPING_IND = "N"
	REC_DQ1_DNA_TYPING_IND = "N"
	REC_DQ2_DNA_TYPING_IND = "N"
	REC_DP1_DNA_TYPING_IND = "N"
	REC_DP2_DNA_TYPING_IND = "N"
	if (re.search(":",REC_A1)):
		REC_A1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_A2)):
		REC_A2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_C1)):
		REC_C1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_C2)):
		REC_C2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_B1)):
		REC_B1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_B2)):
		REC_B2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DR1)):
		REC_DR1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DR2)):
		REC_DR2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DQ1)):
		REC_DQ1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DQ2)):
		REC_DQ2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DP1)):
		REC_DP1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DP2)):
		REC_DP2_DNA_TYPING_IND = "Y"

	# skip donor-recipient pair if either donor or recipient is missing A or B or DR
	if (((DON_A1 == "") & (DON_A2 == "")) | ((DON_B1 == "") & (DON_B2 == "")) | ((DON_DR1 == "") & (DON_DR2 == "")) | 
		((REC_A1 == "") & (REC_A2 == "")) | ((REC_B1 == "") & (REC_B2 == "")) | ((REC_DR1 == "") & (REC_DR2 == ""))):
		if ((DON_A1 == "") & (DON_A2 == "")):
			print ("Missing A in donor: " + str(getattr(row, 'PX_ID')))
		if ((DON_B1 == "") & (DON_B2 == "")):
			print ("Missing B in donor: " + str(getattr(row, 'PX_ID')))
		if ((DON_DR1 == "") & (DON_DR2 == "")):
			print ("Missing DR in donor: " + str(getattr(row, 'PX_ID')))
		if ((REC_A1 == "") & (REC_A2 == "")):
			print ("Missing A in recip: " + str(getattr(row, 'PX_ID')))
		if ((REC_B1 == "") & (REC_B2 == "")):
			print ("Missing B in recip: " + str(getattr(row, 'PX_ID')))
		if ((REC_DR1 == "") & (REC_DR2 == "")):
				print ("Missing DR in recip: " + str(getattr(row, 'PX_ID')))
		# print ("Missing A or B or DR in donor or recip: " + str(getattr(row, 'PX_ID')))
		print (row)
		continue

	# expand_antigens()


	# print GLID file
	# print PULL file

	# PX_ID is donor_ID
	PX_ID = str(getattr(row, 'PX_ID'))

	# print ("A " + DON_A1 + " " + DON_A2 + " DNA " + DON_A1_DNA_TYPING_IND + " " + DON_A2_DNA_TYPING_IND)

	don_gl_A = loc_gl("A",DON_A1,DON_A2,DON_A1_DNA_TYPING_IND,DON_A2_DNA_TYPING_IND,XX_hash)
	don_gl_B = loc_gl("B",DON_B1,DON_B2,DON_B1_DNA_TYPING_IND,DON_B2_DNA_TYPING_IND,XX_hash)
	don_gl_C = loc_gl("C",DON_C1,DON_C2,DON_C1_DNA_TYPING_IND,DON_C2_DNA_TYPING_IND,XX_hash)
	don_gl_DRB1 = loc_gl("DRB1",DON_DR1,DON_DR2,DON_DR1_DNA_TYPING_IND,DON_DR2_DNA_TYPING_IND,XX_hash)
	don_gl_DQB1 = loc_gl("DQB1",DON_DQ1,DON_DQ2,DON_DQ1_DNA_TYPING_IND,DON_DQ2_DNA_TYPING_IND,XX_hash)
	don_gl_DPB1 = loc_gl("DPB1",DON_DP1,DON_DP2,DON_DP1_DNA_TYPING_IND,DON_DP2_DNA_TYPING_IND,XX_hash)
	# print (PX_ID, don_gl_A, don_gl_B, don_gl_C, don_gl_DRB1, don_gl_DQB1, don_gl_DPB1)

	rec_gl_A = loc_gl("A",REC_A1,REC_A2,REC_A1_DNA_TYPING_IND,REC_A2_DNA_TYPING_IND,XX_hash)
	rec_gl_B = loc_gl("B",REC_B1,REC_B2,REC_B1_DNA_TYPING_IND,REC_B2_DNA_TYPING_IND,XX_hash)
	rec_gl_C = loc_gl("C",REC_C1,REC_C2,REC_C1_DNA_TYPING_IND,REC_C2_DNA_TYPING_IND,XX_hash)
	rec_gl_DRB1 = loc_gl("DRB1",REC_DR1,REC_DR2,REC_DR1_DNA_TYPING_IND,REC_DR2_DNA_TYPING_IND,XX_hash)
	rec_gl_DQB1 = loc_gl("DQB1",REC_DQ1,REC_DQ2,REC_DQ1_DNA_TYPING_IND,REC_DQ2_DNA_TYPING_IND,XX_hash)
	rec_gl_DPB1 = loc_gl("DPB1",REC_DP1,REC_DP2,REC_DP1_DNA_TYPING_IND,REC_DP2_DNA_TYPING_IND,XX_hash)
	# print (PX_ID, rec_gl_A, rec_gl_B, rec_gl_C, rec_gl_DRB1, rec_gl_DQB1, rec_gl_DPB1)

	DON_ID = "D" + PX_ID
	REC_ID = "R" + PX_ID

	don_glid_A = 0
	don_glid_B = 0
	don_glid_C = 0
	don_glid_DRBX = 0
	don_glid_DRB1 = 0
	don_glid_DQA1 = 0
	don_glid_DQB1 = 0
	don_glid_DPA1 = 0
	don_glid_DPB1 = 0
	rec_glid_A = 0
	rec_glid_B = 0
	rec_glid_C = 0
	rec_glid_DRBX = 0
	rec_glid_DRB1 = 0
	rec_glid_DQA1 = 0
	rec_glid_DQB1 = 0
	rec_glid_DPA1 = 0
	rec_glid_DPB1 = 0

	if (don_gl_A != "SKIP"):
		if don_gl_A not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_A] = glid_index
			don_glid_A = glid_index
		else:
			don_glid_A = glid_dict[don_gl_A]
	else:
		print ("Invalid A in donor: " + str(getattr(row, 'PX_ID')) + " A " + DON_A1 + " " + DON_A2)
		print (row)
		continue
	if (don_gl_B != "SKIP"):
		if don_gl_B not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_B] = glid_index
			don_glid_B = glid_index
		else:
			don_glid_B = glid_dict[don_gl_B]
	else:
		print ("Invalid B in donor: " + str(getattr(row, 'PX_ID')) + " B " + DON_B1 + " " + DON_B2)
		print (row)
		continue 
	if (don_gl_C != "SKIP"):
		if don_gl_C not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_C] = glid_index
			don_glid_C = glid_index
		else:
			don_glid_C = glid_dict[don_gl_C]
	if (don_gl_DRB1 != "SKIP"):
		if don_gl_DRB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DRB1] = glid_index
			don_glid_DRB1 = glid_index
		else:
			don_glid_DRB1 = glid_dict[don_gl_DRB1]
	else:
		print ("Invalid DR in donor: " + str(getattr(row, 'PX_ID')) + " DR " + DON_DR1 + " " + DON_DR2)
		print (row)
		continue
	if (don_gl_DQB1 != "SKIP"):
		if don_gl_DQB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DQB1] = glid_index
			don_glid_DQB1 = glid_index
		else:
			don_glid_DQB1 = glid_dict[don_gl_DQB1]
	if (don_gl_DPB1 != "SKIP"):
		if don_gl_DPB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DPB1] = glid_index
			don_glid_DPB1 = glid_index
		else:
			don_glid_DPB1 = glid_dict[don_gl_DPB1]

	if (rec_gl_A != "SKIP"):
		if rec_gl_A not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_A] = glid_index
			rec_glid_A = glid_index
		else:
			rec_glid_A = glid_dict[rec_gl_A]
	else:
		print ("Invalid A in recip: " + str(getattr(row, 'PX_ID')) + " A " + REC_A1 + " " + REC_A2)
		print (row)
		continue
	if (rec_gl_B != "SKIP"):
		if rec_gl_B not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_B] = glid_index
			rec_glid_B = glid_index
		else:
			rec_glid_B = glid_dict[rec_gl_B]
	else:
		print ("Invalid B in recip: " + str(getattr(row, 'PX_ID')) + " B " + REC_B1 + " " + REC_B2)
		print (row)
		continue
	if (rec_gl_C != "SKIP"):
		if rec_gl_C not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_C] = glid_index
			rec_glid_C = glid_index
		else:
			rec_glid_C = glid_dict[rec_gl_C]
	if (rec_gl_DRB1 != "SKIP"):
		if rec_gl_DRB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DRB1] = glid_index
			rec_glid_DRB1 = glid_index
		else:
			rec_glid_DRB1 = glid_dict[rec_gl_DRB1]
	else:
		print ("Invalid DR in recip: " + str(getattr(row, 'PX_ID')) + " DR " + REC_DR1 + " " + REC_DR2)
		print (row)
		continue
	if (rec_gl_DQB1 != "SKIP"):
		if rec_gl_DQB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DQB1] = glid_index
			rec_glid_DQB1 = glid_index
		else:
			rec_glid_DQB1 = glid_dict[rec_gl_DQB1]
	if (rec_gl_DPB1 != "SKIP"):
		if rec_gl_DPB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DPB1] = glid_index
			rec_glid_DPB1 = glid_index
		else:
			rec_glid_DPB1 = glid_dict[rec_gl_DPB1]


	pull_string = DON_ID + "," + DON_RACE + "," + ','.join([str(don_glid_A),str(don_glid_B),str(don_glid_C),str(don_glid_DRB1),str(don_glid_DRBX),str(don_glid_DQA1),str(don_glid_DQB1),str(don_glid_DPA1),str(don_glid_DPB1)])
	PULL.write(pull_string + "\n")
	pull_string = REC_ID + "," + CAN_RACE + "," + ','.join([str(rec_glid_A),str(rec_glid_B),str(rec_glid_C),str(rec_glid_DRB1),str(rec_glid_DRBX),str(rec_glid_DQA1),str(rec_glid_DQB1),str(rec_glid_DPA1),str(rec_glid_DPB1)])
	PULL.write(pull_string + "\n")

PULL.close()

# print GLID file
for gl in glid_dict:
	GLID.write(str(glid_dict[gl]) + "," + gl + "\n")

GLID.close()

	