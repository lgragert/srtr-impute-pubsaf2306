#!/usr/bin/env python
import pandas
import numpy
import re
import requests
# import pyard
import sys

# initialize pyARD and roll up typing to 'g' group for imputation
# XX codes are always in latest IMGT/HLA allele list
# ard = pyard.init()

# global XX code GL hash
XX_hash = {}
dna2ser = {}
seroGL = {}

sysarg_donor = sys.argv[1]
sysarg_tx_ty = sys.argv[2]
print("DON Type:", sysarg_donor)
print("TX Type:", sysarg_tx_ty)

donor_type = str(sysarg_donor)  # deceased or living
tx_type = str(sysarg_tx_ty)     # heart_lungs, heart, lungs, kidney_pancreas, kidney, pancreas, liver, intestines

# donor_type = "living"
# tx_type = "kidney"

def loc_gl (loc, typ1, typ2, dna1, dna2, XX_hash):

	# manage homozygotes by copying over antigen
	# TODO - modify to add null allele list when homozygous
	if (typ2 == ""): 
		typ2 = typ1
		dna2 = dna1
	if (typ1 == ""):
		typ1 = typ2
		dna1 = dna2
	if ((typ1 == "") | (typ2 == "")):
		if (loc == "DRB345"):
			return '+'.join([null_alleles_glstring_locus["DRB345"],null_alleles_glstring_locus["DRB345"]])
		# if (loc in ["A","B","DRB1"]):
		# 	print ("Either no typ1 or typ2 in loc_gl: " + str(getattr(row, 'PX_ID')) + " Typ1: " + typ1 + " Typ2: " + typ2 + "^")
		return "SKIP"

	# print (loc + " TYP: " + typ1 + " " + typ2 + " DNA: " + dna)	
	
	al1 = ""
	al2 = ""

	# TODO - handle C13
	# if ((typ1 == "13:XX") | (typ2 == "13:XX")):
		# Cw*1301	Sequence shown to be in error (July 2002) - https://hla.alleles.org/alleles/deleted.html
		# print ("Deleted allele C*13:XX in HLA typing: " + str(getattr(row, 'PX_ID')))
	# 	return "SKIP"

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
		if (loc == "DRB345"):
			loc_ser = "DR"
			if ("*" in typ1):
				loc_ser = "DRB"		
		if (loc == "DQA1"):
			loc_ser = "DQA1*"
		if (loc == "DQB1"):
			loc_ser = "DQ"
		if (loc == "DPA1"):
			loc_ser = "DPA1*"
		if (loc == "DPB1"):
			loc_ser = "DP"
			# print ("Legacy DP antigen-level typing  - not using for imputation:" + str(getattr(row, 'PX_ID')))
			return "SKIP"  # DP serology is inactive
		loctyp1 = loc_ser + typ1
		# print ("SERO " + loctyp1 + "" + loctyp2)
		al1 = seroGL[loctyp1]

	# dna typing
	elif (dna1 == "Y"):
		loc_dna = loc
		typ1_fields = typ1.split(':')
		if (loc == "DRB345"):
			if ("*" in typ1): # DRB345
				loc_dna = "DRB"	
				loctyp1 = loc_dna + typ2
		else:
			subtype1 = typ1_fields[1]
			loctyp1 = '*'.join([loc_dna,typ1])
		al1 = loctyp1
		# two-field analysis - ARD rollup no longer being used
		# al1 = ard.redux(al1,'lgx')

	if (dna2 == "N"):
		loc_ser = loc
		if (loc == "DRB1"):
			loc_ser = "DR"
		if (loc == "DRB345"):
			loc_ser = "DR"
			if ("*" in typ2):
				loc_ser = "DRB"
		if (loc == "DQA1"):
			loc_ser = "DQA1*"
		if (loc == "DQB1"):
			loc_ser = "DQ"
		if (loc == "DPA1"):
			loc_ser = "DPA1*"
		if (loc == "DPB1"):
			loc_ser = "DP"
			# print ("Legacy DP antigen-level typing - not using for imputation: " + str(getattr(row, 'PX_ID')))
			return "SKIP"  # DP serology is inactive
		loctyp2 = loc_ser + typ2
		# 2nd typing can indicate either homozygous or null
		if (typ1 == typ2): 
			al2 = seroGL[loctyp1] + '/' + null_alleles_glstring_locus[loc]
		else:
			al2 = seroGL[loctyp2]

	# dna typing
	elif (dna2 == "Y"):
		loc_dna = loc
		if (typ2 == ""): 
			al2 = loctyp1 + null_alleles_glstring_locus[loc]
		else:
			typ2_fields = typ2.split(':')
			if (loc_dna == "DRB345"):
				if ("*" in typ2): # DRB345
					loc_dna = "DRB"
					loctyp2 = loc_dna + typ2
			else:
				subtype2 = typ2_fields[1]
				loctyp2 = '*'.join([loc_dna,typ2])
			al2 = loctyp2
		# # two-field analysis - ARD rollup no longer being used
		# al2 = ard.redux(al2,'lgx')

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

unos_sero_filename = "OPTN_antigens_to_alleles_CPRA.txt"
unos_sero_df = pandas.read_csv(unos_sero_filename, sep='\t')


# get list of two-field alleles that include nulls from WMDA file
wmda_filename  = "rel_dna_ser.txt"
column_names = ['Who_Locus', 'Who_Allele_Name', 'Who_Unambiguous_Ag', 'Possible_Ag', 'Assumed_Ag', 'Expert_Assigned_Ag']
wmda_df = pandas.read_csv(wmda_filename, sep=';', skiprows=6, header=None, names=column_names)

wmda_file = open(wmda_filename, 'r')
null_alleles_locus = {} # list of null alleles per locus
for line in wmda_file:
	# skip header and some nonclassical HLA loci
	if ((line.startswith('#')) or (line.startswith('DMA')) or (line.startswith('DMB')) or (line.startswith('DOA')) or
		(line.startswith('DOB')) or (line.startswith('DPA2')) or (line.startswith('DQA2')) or (line.startswith('DPB')) or
		(line.startswith('DRA')) or	(line.startswith('DRB2'))):
		continue

	# skip rest of file after DRB6
	if (line.startswith('DRB6')):
		break

	line = line.strip('\n')
	(who_locus,who_allele_name,who_unambiguous_ag,possible_ag,assumed_ag,expert_assigned_ag) = line.split(';')

	imgt_hla_allele = who_locus + who_allele_name
	allele_fields = who_allele_name.split(':')
	allele_group = allele_fields[0]
	protein = allele_fields[1]
	has_expression_char = re.search('[A-Z]$', who_allele_name)
	imgt_two_field_allele = who_locus + allele_group + ":" + protein # no shortnulls in freqs

	# trim off trailing "*" character from locus name
	who_locus = who_locus[:-1]
	if who_locus in ["DRB3","DRB4","DRB5"]:
		who_locus = "DRB345"

	expression_char = ""
	if (has_expression_char):
		expression_char = has_expression_char.group(0)
	if (expression_char == "N"): # add two-field null allele to locus-specific list
		if who_locus in null_alleles_locus.keys():
			if imgt_two_field_allele in null_alleles_locus[who_locus]:
				continue
			else:
				null_alleles_locus[who_locus].append(imgt_two_field_allele)
		else:
			null_alleles_locus[who_locus] = [imgt_two_field_allele]

null_alleles_locus["DRB345"].append("DRBX*NNNN")

null_alleles_glstring_locus = {}

for locus in null_alleles_locus:
	null_alleles_glstring_locus[locus] = '/'.join(null_alleles_locus[locus])
	print(null_alleles_glstring_locus[locus])


# load file generated from HLA merge
# tx_ki_hla_filename = "tx_ki_hla_9loc.csv"
# tx_ki_hla_filename = "deceased_kidney_hla_9loc.csv"
tx_ki_hla_filename = f"{donor_type}_{tx_type}_hla_9loc.csv"
tx_ki_hla = pandas.read_csv(tx_ki_hla_filename, index_col=None, encoding='Latin-1', low_memory=False, dtype=str)

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
	# if (CAN_RACE == "1024: Unknown (for Donor Referral only)"):
	# 	CAN_RACE = "UNK"
	# if (DON_RACE == "1024: Unknown (for Donor Referral only)"):
	# 	DON_RACE = "UNK"	
	# if (CAN_RACE == "Missing"):
	# 	CAN_RACE = "UNK"
	# if (DON_RACE == "Missing"):
	#	DON_RACE = "UNK"	

	# handle nan
	DON_A1 = str(getattr(row,'DON_A1'))
	DON_A2 = str(getattr(row,'DON_A2'))
	DON_C1 = str(getattr(row,'DON_C1'))
	DON_C2 = str(getattr(row,'DON_C2'))
	DON_B1 = str(getattr(row,'DON_B1'))
	DON_B2 = str(getattr(row,'DON_B2'))
	DON_DR1 = str(getattr(row,'DON_DR1'))
	DON_DR2 = str(getattr(row,'DON_DR2'))
	DON_DQA1 = str(getattr(row,'DON_DQA1'))
	DON_DQA2 = str(getattr(row,'DON_DQA2'))
	DON_DQ1 = str(getattr(row,'DON_DQ1'))
	DON_DQ2 = str(getattr(row,'DON_DQ2'))
	DON_DPA1 = str(getattr(row,'DON_DPA1'))
	DON_DPA2 = str(getattr(row,'DON_DPA2'))
	DON_DP1 = str(getattr(row,'DON_DP1'))
	DON_DP2 = str(getattr(row,'DON_DP2'))
	# DON_DR51 = str(getattr(row,'DON_DR51')) # this field is only Negative or Positive
	# DON_DR52 = str(getattr(row,'DON_DR52')) # this field is only Negative or Positive
	# DON_DR53 = str(getattr(row,'DON_DR53')) # this field is only Negative or Positive
	DON_DRB3_1 = str(getattr(row,'DON_DRB3_1'))
	DON_DRB3_2 = str(getattr(row,'DON_DRB3_2'))
	DON_DRB4_1 = str(getattr(row,'DON_DRB4_1'))
	DON_DRB4_2 = str(getattr(row,'DON_DRB4_2'))
	DON_DRB5_1 = str(getattr(row,'DON_DRB5_1'))
	DON_DRB5_2 = str(getattr(row,'DON_DRB5_2'))

	# 97	97: Unknown
	# 98	98: No second antigen detected
	# 99	99: Not Tested

	if ((DON_A1 == "nan") or (DON_A1 == "Null_allele") or (DON_A1 == "97") or (DON_A1 == "98") or (DON_A1 == "99")):
		DON_A1 = ""
	if (DON_A2 == "nan") or (DON_A2 == "Null_allele") or (DON_A2 == "97") or (DON_A2 == "98") or (DON_A2 == "99"):
		DON_A2 = ""
	if (DON_C1 == "nan") or (DON_C1 == "Null_allele") or (DON_C1 == "97") or (DON_C1 == "98") or (DON_C1 == "99") or (DON_C1 == "No Ag detected"):
		DON_C1 = ""
	if (DON_C2 == "nan") or (DON_C2 == "Null_allele") or (DON_C2 == "97") or (DON_C2 == "98") or (DON_C2 == "99") or (DON_C2 == "No Ag detected"):
		DON_C2 = ""
	if (DON_B1 == "nan") or (DON_B1 == "Null_allele") or (DON_B1 == "97") or (DON_B1 == "98") or (DON_B1 == "99"):
		DON_B1 = ""
	if (DON_B2 == "nan") or (DON_B2 == "Null_allele") or (DON_B2 == "97") or (DON_B2 == "98") or (DON_B2 == "99"):
		DON_B2 = ""
	if (DON_DR1 == "nan") or (DON_DR1 == "Null_allele") or (DON_DR1 == "97") or (DON_DR1 == "98") or (DON_DR1 == "99"):
		DON_DR1 = ""
	if (DON_DR2 == "nan") or (DON_DR2 == "Null_allele") or (DON_DR2 == "97") or (DON_DR2 == "98") or (DON_DR2 == "99"):
		DON_DR2 = ""
	if (DON_DQA1 == "nan") or (DON_DQA1 == "Null_allele") or (DON_DQA1 == "97") or (DON_DQA1 == "98") or (DON_DQA1 == "99"):
		DON_DQA1 = ""
	if (DON_DQA2 == "nan") or (DON_DQA2 == "Null_allele") or (DON_DQA2 == "97") or (DON_DQA2 == "98") or (DON_DQA2 == "99"):
		DON_DQA2 = ""
	if (DON_DQ1 == "nan") or (DON_DQ1 == "Null_allele") or (DON_DQ1 == "97") or (DON_DQ1 == "98") or (DON_DQ1 == "99"):
		DON_DQ1 = ""
	if (DON_DQ2 == "nan") or (DON_DQ2 == "Null_allele") or (DON_DQ2 == "97") or (DON_DQ2 == "98") or (DON_DQ2 == "99"):
		DON_DQ2 = ""
	if (DON_DPA1 == "nan") or (DON_DPA1 == "Null_allele") or (DON_DPA1 == "97") or (DON_DPA1 == "98") or (DON_DPA1 == "99"):
		DON_DPA1 = ""
	if (DON_DPA2 == "nan") or (DON_DPA2 == "Null_allele") or (DON_DPA2 == "97") or (DON_DPA2 == "98") or (DON_DPA2 == "99"):
		DON_DPA2 = ""
	if (DON_DP1 == "nan") or (DON_DP1 == "Null_allele") or (DON_DP1 == "97") or (DON_DP1 == "98") or (DON_DP1 == "99"):
		DON_DP1 = ""
	if (DON_DP2 == "nan") or (DON_DP2 == "Null_allele") or (DON_DP2 == "97") or (DON_DP2 == "98") or (DON_DP2 == "99"):
		DON_DP2 = ""
	if (DON_DRB3_1 == "nan") or (DON_DRB3_1 == "Null_allele"):
		DON_DRB3_1 = ""
	if (DON_DRB3_2 == "nan") or (DON_DRB3_2 == "Null_allele"):
		DON_DRB3_2 = ""
	if (DON_DRB4_1 == "nan") or (DON_DRB4_1 == "Null_allele"):
		DON_DRB4_1 = ""
	if (DON_DRB4_2 == "nan") or (DON_DRB4_2 == "Null_allele"):
		DON_DRB4_2 = ""
	if (DON_DRB5_1 == "nan") or (DON_DRB5_1 == "Null_allele"):
		DON_DRB5_1 = ""
	if (DON_DRB5_2 == "nan") or (DON_DRB5_2 == "Null_allele"):
		DON_DRB5_2 = ""


	DON_A1_DNA_TYPING_IND = "N"
	DON_A2_DNA_TYPING_IND = "N"
	DON_C1_DNA_TYPING_IND = "N"
	DON_C2_DNA_TYPING_IND = "N"
	DON_B1_DNA_TYPING_IND = "N"
	DON_B2_DNA_TYPING_IND = "N"
	DON_DR1_DNA_TYPING_IND = "N"
	DON_DR2_DNA_TYPING_IND = "N"
	DON_DQA1_DNA_TYPING_IND = "N"
	DON_DQA2_DNA_TYPING_IND = "N"
	DON_DQ1_DNA_TYPING_IND = "N"
	DON_DQ2_DNA_TYPING_IND = "N"
	DON_DPA1_DNA_TYPING_IND = "N"
	DON_DPA2_DNA_TYPING_IND = "N"
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
	if (re.search(":",DON_DQA1)):
		DON_DQA1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DQA2)):
		DON_DQA2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DQ1)):
		DON_DQ1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DQ2)):
		DON_DQ2_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DPA1)):
		DON_DPA1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DPA2)):
		DON_DPA2_DNA_TYPING_IND = "Y"
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
	REC_DQA1 = str(getattr(row,'REC_DQA1'))
	REC_DQA2 = str(getattr(row,'REC_DQA2'))
	REC_DQ1 = str(getattr(row,'REC_DQW1'))
	REC_DQ2 = str(getattr(row,'REC_DQW2'))
	REC_DPA1 = str(getattr(row,'REC_DPA1'))
	REC_DPA2 = str(getattr(row,'REC_DPA2'))
	REC_DP1 = str(getattr(row,'REC_DPW1'))
	REC_DP2 = str(getattr(row,'REC_DPW2'))
	# REC_DR51 = str(getattr(row,'REC_DRW51')) # this field is only Negative or Positive
	# REC_DR52 = str(getattr(row,'REC_DRW52')) # this field is only Negative or Positive
	# REC_DR53 = str(getattr(row,'REC_DRW53')) # this field is only Negative or Positive
	REC_DRB3_1 = str(getattr(row,'REC_DRB3_1'))
	REC_DRB3_2 = str(getattr(row,'REC_DRB3_2'))
	REC_DRB4_1 = str(getattr(row,'REC_DRB4_1'))
	REC_DRB4_2 = str(getattr(row,'REC_DRB4_2'))
	REC_DRB5_1 = str(getattr(row,'REC_DRB5_1'))
	REC_DRB5_2 = str(getattr(row,'REC_DRB5_2'))

	if ((REC_A1 == "nan") or (REC_A1 == "Null_allele") or (REC_A1 == "97") or (REC_A1 == "98") or (REC_A1 == "99")):
		REC_A1 = ""
	if (REC_A2 == "nan") or (REC_A2 == "Null_allele") or (REC_A2 == "97") or (REC_A2 == "98") or (REC_A2 == "99"):
		REC_A2 = ""
	if (REC_C1 == "nan") or (REC_C1 == "Null_allele") or (REC_C1 == "97") or (REC_C1 == "98") or (REC_C1 == "99") or (REC_C1 == "No Ag detected"):
		REC_C1 = ""
	if (REC_C2 == "nan") or (REC_C2 == "Null_allele") or (REC_C2 == "97") or (REC_C2 == "98") or (REC_C2 == "99") or (REC_C2 == "No Ag detected"):
		REC_C2 = ""
	if (REC_B1 == "nan") or (REC_B1 == "Null_allele") or (REC_B1 == "97") or (REC_B1 == "98") or (REC_B1 == "99"):
		REC_B1 = ""
	if (REC_B2 == "nan") or (REC_B2 == "Null_allele") or (REC_B2 == "97") or (REC_B2 == "98") or (REC_B2 == "99"):
		REC_B2 = ""
	if (REC_DR1 == "nan") or (REC_DR1 == "Null_allele") or (REC_DR1 == "97") or (REC_DR1 == "98") or (REC_DR1 == "99"):
		REC_DR1 = ""
	if (REC_DR2 == "nan") or (REC_DR2 == "Null_allele") or (REC_DR2 == "97") or (REC_DR2 == "98") or (REC_DR2 == "99"):
		REC_DR2 = ""
	if (REC_DQA1 == "nan") or (REC_DQA1 == "Null_allele") or (REC_DQA1 == "97") or (REC_DQA1 == "98") or (REC_DQA1 == "99"):
		REC_DQA1 = ""
	if (REC_DQA2 == "nan") or (REC_DQA2 == "Null_allele") or (REC_DQA2 == "97") or (REC_DQA2 == "98") or (REC_DQA2 == "99"):
		REC_DQA2 = ""
	if (REC_DQ1 == "nan") or (REC_DQ1 == "Null_allele") or (REC_DQ1 == "97") or (REC_DQ1 == "98") or (REC_DQ1 == "99") or (len(REC_DQ1) == 2): # TODO - 2-digit typing not found in seroGL[DQ]
		REC_DQ1 = ""
	if (REC_DQ2 == "nan") or (REC_DQ2 == "Null_allele") or (REC_DQ2 == "97") or (REC_DQ2 == "98") or (REC_DQ2 == "99") or (len(REC_DQ2) == 2): # TODO - 2-digit typing not found in seroGL[DQ]
		REC_DQ2 = ""
	if (REC_DPA1 == "nan") or (REC_DPA1 == "Null_allele") or (REC_DPA1 == "97") or (REC_DPA1 == "98") or (REC_DPA1 == "99"):
		REC_DPA1 = ""
	if (REC_DPA2 == "nan") or (REC_DPA2 == "Null_allele") or (REC_DPA2 == "97") or (REC_DPA2 == "98") or (REC_DPA2 == "99"):
		REC_DPA2 = ""
	if (REC_DP1 == "nan") or (REC_DP1 == "Null_allele") or (REC_DP1 == "97") or (REC_DP1 == "98") or (REC_DP1 == "99"):
		REC_DP1 = ""
	if (REC_DP2 == "nan") or (REC_DP2 == "Null_allele") or (REC_DP2 == "97") or (REC_DP2 == "98") or (REC_DP2 == "99"):
		REC_DP2 = ""
	if (REC_DRB3_1 == "nan") or (REC_DRB3_1 == "Null_allele"):
		REC_DRB3_1 = ""
	if (REC_DRB3_2 == "nan") or (REC_DRB3_2 == "Null_allele"):
		REC_DRB3_2 = ""
	if (REC_DRB4_1 == "nan") or (REC_DRB4_1 == "Null_allele"):
		REC_DRB4_1 = ""
	if (REC_DRB4_2 == "nan") or (REC_DRB4_2 == "Null_allele"):
		REC_DRB4_2 = ""
	if (REC_DRB5_1 == "nan") or (REC_DRB5_1 == "Null_allele"):
		REC_DRB5_1 = ""
	if (REC_DRB5_2 == "nan") or (REC_DRB5_2 == "Null_allele"):
		REC_DRB5_2 = ""

	REC_A1_DNA_TYPING_IND = "N"
	REC_A2_DNA_TYPING_IND = "N"
	REC_C1_DNA_TYPING_IND = "N"
	REC_C2_DNA_TYPING_IND = "N"
	REC_B1_DNA_TYPING_IND = "N"
	REC_B2_DNA_TYPING_IND = "N"
	REC_DR1_DNA_TYPING_IND = "N"
	REC_DR2_DNA_TYPING_IND = "N"
	REC_DQA1_DNA_TYPING_IND = "N"
	REC_DQA2_DNA_TYPING_IND = "N"
	REC_DQ1_DNA_TYPING_IND = "N"
	REC_DQ2_DNA_TYPING_IND = "N"
	REC_DPA1_DNA_TYPING_IND = "N"
	REC_DPA2_DNA_TYPING_IND = "N"
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
	if (re.search(":",REC_DQA1)):
		REC_DQA1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DQA2)):
		REC_DQA2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DQ1)):
		REC_DQ1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DQ2)):
		REC_DQ2_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DPA1)):
		REC_DPA1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DPA2)):
		REC_DPA2_DNA_TYPING_IND = "Y"
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


	# handle DRB3/4/5
	# Previously REC_DR51, REC_DR52, REC_DR53, DON_DR51, DON_DR52, DON_DR52
	DON_DRB345_1 = ""
	DON_DRB345_2 = ""	
	if (DON_DRB3_1 != ""):
		if (DON_DRB345_1 == ""):
			DON_DRB345_1 = DON_DRB3_1
		else:
			DON_DRB345_2 = DON_DRB3_1
	if (DON_DRB3_2 != ""):
		if (DON_DRB345_1 == ""):
			DON_DRB345_1 = DON_DRB3_2
		else:
			DON_DRB345_2 = DON_DRB3_2
	if (DON_DRB4_1 != ""):
		if (DON_DRB345_1 == ""):
			DON_DRB345_1 = DON_DRB4_1
		else:
			DON_DRB345_2 = DON_DRB4_1
	if (DON_DRB4_2 != ""):
		if (DON_DRB345_1 == ""):
			DON_DRB345_1 = DON_DRB4_2
		else:
			DON_DRB345_2 = DON_DRB4_2
	if (DON_DRB5_1 != ""):
		if (DON_DRB345_1 == ""):
			DON_DRB345_1 = DON_DRB5_1
		else:
			DON_DRB345_2 = DON_DRB5_1
	if (DON_DRB5_2 != ""):
		if (DON_DRB345_1 == ""):
			DON_DRB345_1 = DON_DRB5_2
		else:
			DON_DRB345_2 = DON_DRB5_2

	DON_DRB345_1_DNA_TYPING_IND = "N"
	DON_DRB345_2_DNA_TYPING_IND = "N"
	if (re.search(":",DON_DRB345_1)):
		DON_DRB345_1_DNA_TYPING_IND = "Y"
	if (re.search(":",DON_DRB345_2)):
		DON_DRB345_2_DNA_TYPING_IND = "Y"


	REC_DRB345_1 = ""
	REC_DRB345_2 = ""	
	if (REC_DRB3_1 != ""):
		if (REC_DRB345_1 == ""):
			REC_DRB345_1 = REC_DRB3_1
		else:
			REC_DRB345_2 = REC_DRB3_1
	if (REC_DRB3_2 != ""):
		if (REC_DRB345_1 == ""):
			REC_DRB345_1 = REC_DRB3_2
		else:
			REC_DRB345_2 = REC_DRB3_2
	if (REC_DRB4_1 != ""):
		if (REC_DRB345_1 == ""):
			REC_DRB345_1 = REC_DRB4_1
		else:
			REC_DRB345_2 = REC_DRB4_1
	if (REC_DRB4_2 != ""):
		if (REC_DRB345_1 == ""):
			REC_DRB345_1 = REC_DRB4_2
		else:
			REC_DRB345_2 = REC_DRB4_2
	if (REC_DRB5_1 != ""):
		if (REC_DRB345_1 == ""):
			REC_DRB345_1 = REC_DRB5_1
		else:
			REC_DRB345_2 = REC_DRB5_1
	if (REC_DRB5_2 != ""):
		if (REC_DRB345_1 == ""):
			REC_DRB345_1 = REC_DRB5_2
		else:
			REC_DRB345_2 = REC_DRB5_2

	REC_DRB345_1_DNA_TYPING_IND = "N"
	REC_DRB345_2_DNA_TYPING_IND = "N"
	if (re.search(":",REC_DRB345_1)):
		REC_DRB345_1_DNA_TYPING_IND = "Y"
	if (re.search(":",REC_DRB345_2)):
		REC_DRB345_2_DNA_TYPING_IND = "Y"

	# KeyError: 'C1'
	# If DON or REC has C1 or C2 that is single-digit, then prepend "0" to the single digit
	if (DON_C1 != ""):
		if (len(DON_C1) == 1):
			DON_C1 = "0" + DON_C1
	if (DON_C2 != ""):
		if (len(DON_C2) == 1):
			DON_C2 = "0" + DON_C2
	if (REC_C1 != ""):
		if (len(REC_C1) == 1):
			REC_C1 = "0" + REC_C1
	if (REC_C2 != ""):
		if (len(REC_C2) == 1):
			REC_C2 = "0" + REC_C2

	# KeyError: 'DRHLA-DR103'
	# If DON or REC has DR1 or DR2 that is 'HLA-DR103', then change it to '103'
	if (DON_DR1 != ""):
		if (DON_DR1 == "HLA-DR103"):
			DON_DR1 = "103"
	if (DON_DR2 != ""):
		if (DON_DR2 == "HLA-DR103"):
			DON_DR2 = "103"
	if (REC_DR1 != ""):
		if (REC_DR1 == "HLA-DR103"):
			REC_DR1 = "103"
	if (REC_DR2 != ""):
		if (REC_DR2 == "HLA-DR103"):
			REC_DR2 = "103"
	
	# REC DQW 2-digit typing not found in seroGL[DQ]
	# If REC_HISTO_TX_ID is found in ['1045249', '1058933', '1075858', '1078954', '1080256'], and if REC_DQ1 or REC_DQ2 is 2-digit typing, then change it to ""
	REC_HISTO_TX_ID = str(getattr(row,'REC_HISTO_TX_ID'))
	if REC_HISTO_TX_ID in ['1045249', '1058933', '1075858', '1078954', '1080256']:
		if (len(REC_DQ1) == 2) or (len(REC_DQ2) == 2):
			# Change REC_DQ1 and REC_DQ2 to ""
			REC_DQ1 = ""
			REC_DQ2 = ""

	# print GLID file
	# print PULL file

	# PX_ID is donor_ID
	PX_ID = str(getattr(row, 'PX_ID'))
	# print (PX_ID)
	# print (row)

	# print ("A " + DON_A1 + " " + DON_A2 + " DNA " + DON_A1_DNA_TYPING_IND + " " + DON_A2_DNA_TYPING_IND)
	# print ("C " + DON_C1 + " " + DON_C2 + " DNA " + DON_C1_DNA_TYPING_IND + " " + DON_C2_DNA_TYPING_IND)

	don_gl_A = loc_gl("A",DON_A1,DON_A2,DON_A1_DNA_TYPING_IND,DON_A2_DNA_TYPING_IND,XX_hash)
	don_gl_B = loc_gl("B",DON_B1,DON_B2,DON_B1_DNA_TYPING_IND,DON_B2_DNA_TYPING_IND,XX_hash)
	don_gl_C = loc_gl("C",DON_C1,DON_C2,DON_C1_DNA_TYPING_IND,DON_C2_DNA_TYPING_IND,XX_hash)
	don_gl_DRB345 = loc_gl("DRB345",DON_DRB345_1,DON_DRB345_2,DON_DRB345_1_DNA_TYPING_IND,DON_DRB345_2_DNA_TYPING_IND,XX_hash)
	don_gl_DRB1 = loc_gl("DRB1",DON_DR1,DON_DR2,DON_DR1_DNA_TYPING_IND,DON_DR2_DNA_TYPING_IND,XX_hash)
	don_gl_DQA1 = loc_gl("DQA1",DON_DQA1,DON_DQA2,DON_DQA1_DNA_TYPING_IND,DON_DQA2_DNA_TYPING_IND,XX_hash)
	don_gl_DQB1 = loc_gl("DQB1",DON_DQ1,DON_DQ2,DON_DQ1_DNA_TYPING_IND,DON_DQ2_DNA_TYPING_IND,XX_hash)
	don_gl_DPA1 = loc_gl("DPA1",DON_DPA1,DON_DPA2,DON_DPA1_DNA_TYPING_IND,DON_DPA2_DNA_TYPING_IND,XX_hash)
	don_gl_DPB1 = loc_gl("DPB1",DON_DP1,DON_DP2,DON_DP1_DNA_TYPING_IND,DON_DP2_DNA_TYPING_IND,XX_hash)
	# print (PX_ID, don_gl_A, don_gl_B, don_gl_C, don_gl_DRB1, don_gl_DQB1, don_gl_DPB1)

	rec_gl_A = loc_gl("A",REC_A1,REC_A2,REC_A1_DNA_TYPING_IND,REC_A2_DNA_TYPING_IND,XX_hash)
	rec_gl_B = loc_gl("B",REC_B1,REC_B2,REC_B1_DNA_TYPING_IND,REC_B2_DNA_TYPING_IND,XX_hash)
	rec_gl_C = loc_gl("C",REC_C1,REC_C2,REC_C1_DNA_TYPING_IND,REC_C2_DNA_TYPING_IND,XX_hash)
	rec_gl_DRB345 = loc_gl("DRB345",REC_DRB345_1,REC_DRB345_2,REC_DRB345_1_DNA_TYPING_IND,REC_DRB345_2_DNA_TYPING_IND,XX_hash)
	rec_gl_DRB1 = loc_gl("DRB1",REC_DR1,REC_DR2,REC_DR1_DNA_TYPING_IND,REC_DR2_DNA_TYPING_IND,XX_hash)
	rec_gl_DQA1 = loc_gl("DQA1",REC_DQA1,REC_DQA2,REC_DQA1_DNA_TYPING_IND,REC_DQA2_DNA_TYPING_IND,XX_hash)
	rec_gl_DQB1 = loc_gl("DQB1",REC_DQ1,REC_DQ2,REC_DQ1_DNA_TYPING_IND,REC_DQ2_DNA_TYPING_IND,XX_hash)
	rec_gl_DPA1 = loc_gl("DPA1",REC_DPA1,REC_DPA2,REC_DPA1_DNA_TYPING_IND,REC_DPA2_DNA_TYPING_IND,XX_hash)
	rec_gl_DPB1 = loc_gl("DPB1",REC_DP1,REC_DP2,REC_DP1_DNA_TYPING_IND,REC_DP2_DNA_TYPING_IND,XX_hash)
	# print (PX_ID, rec_gl_A, rec_gl_B, rec_gl_C, rec_gl_DRB1, rec_gl_DQB1, rec_gl_DPB1)

	DON_ID = "D" + PX_ID
	REC_ID = "R" + PX_ID

	don_glid_A = 0
	don_glid_B = 0
	don_glid_C = 0
	don_glid_DRB345 = 0
	don_glid_DRB1 = 0
	don_glid_DQA1 = 0
	don_glid_DQB1 = 0
	don_glid_DPA1 = 0
	don_glid_DPB1 = 0
	rec_glid_A = 0
	rec_glid_B = 0
	rec_glid_C = 0
	rec_glid_DRB345 = 0
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
	if (don_gl_DRB345 != "SKIP"):
		if don_gl_DRB345 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DRB345] = glid_index
			don_gl_DRB345 = glid_index
		else:
			don_gl_DRB345 = glid_dict[don_gl_DRB345]
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
	if (don_gl_DQA1 != "SKIP"):
		if don_gl_DQA1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DQA1] = glid_index
			don_glid_DQA1 = glid_index
		else:
			don_glid_DQA1 = glid_dict[don_gl_DQA1]
	if (don_gl_DQB1 != "SKIP"):
		if don_gl_DQB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DQB1] = glid_index
			don_glid_DQB1 = glid_index
		else:
			don_glid_DQB1 = glid_dict[don_gl_DQB1]
	if (don_gl_DPA1 != "SKIP"):
		if don_gl_DPA1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[don_gl_DPA1] = glid_index
			don_glid_DPA1 = glid_index
		else:
			don_glid_DPA1 = glid_dict[don_gl_DPA1]
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
	if (rec_gl_DRB345 != "SKIP"):
		if rec_gl_DRB345 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DRB345] = glid_index
			rec_gl_DRB345 = glid_index
		else:
			rec_gl_DRB345 = glid_dict[rec_gl_DRB345]
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
	if (rec_gl_DQA1 != "SKIP"):
		if rec_gl_DQA1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DQA1] = glid_index
			rec_glid_DQA1 = glid_index
		else:
			rec_glid_DQA1 = glid_dict[rec_gl_DQA1]
	if (rec_gl_DQB1 != "SKIP"):
		if rec_gl_DQB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DQB1] = glid_index
			rec_glid_DQB1 = glid_index
		else:
			rec_glid_DQB1 = glid_dict[rec_gl_DQB1]
	if (rec_gl_DPA1 != "SKIP"):
		if rec_gl_DPA1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DPA1] = glid_index
			rec_glid_DPA1 = glid_index
		else:
			rec_glid_DPA1 = glid_dict[rec_gl_DPA1]
	if (rec_gl_DPB1 != "SKIP"):
		if rec_gl_DPB1 not in glid_dict:
			glid_index = glid_index + 1
			glid_dict[rec_gl_DPB1] = glid_index
			rec_glid_DPB1 = glid_index
		else:
			rec_glid_DPB1 = glid_dict[rec_gl_DPB1]


	pull_string = DON_ID + "," + DON_RACE + "," + ','.join([str(don_glid_A),str(don_glid_B),str(don_glid_C),str(don_glid_DRB1),str(don_glid_DRB345),str(don_glid_DQA1),str(don_glid_DQB1),str(don_glid_DPA1),str(don_glid_DPB1)])
	PULL.write(pull_string + "\n")
	pull_string = REC_ID + "," + CAN_RACE + "," + ','.join([str(rec_glid_A),str(rec_glid_B),str(rec_glid_C),str(rec_glid_DRB1),str(rec_glid_DRB345),str(rec_glid_DQA1),str(rec_glid_DQB1),str(rec_glid_DPA1),str(rec_glid_DPB1)])
	PULL.write(pull_string + "\n")

PULL.close()

# print GLID file
for gl in glid_dict:
	GLID.write(str(glid_dict[gl]) + "," + gl + "\n")

GLID.close()

	
