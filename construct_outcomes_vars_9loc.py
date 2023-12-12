#!/usr/bin/env python3

from os import path
import glob
import os
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime
import time


replicates = 10

for i in range(1, replicates+1):
	# load in raw design matrix
	srtr_AA_MM_filename = "SRTR_AA_MM_9loc_" + str(i) + ".txt"
	design_matrix = pd.read_csv(srtr_AA_MM_filename, sep='\t', index_col=None, encoding='Latin-1', low_memory=False)

	# print (design_matrix)

	# Once the input data are worked with, use post-transplant survival model;

	# Restrict post-transplant data to transplants that fit filter;
	# let ptfilter = Rec_age_at_tx ge 18 and org_ty in ('KI') and REC_TX_ORG_TY = "KI"
	#                 /*and don_ty eq 'C'*/ and don_age ge 9 /*and CAN_SSN_FLG = "V"*/; 
	# let censor_dt = "01JUL2011"D;


	design_matrix['TFL_DEATH_DT'] = pd.to_datetime(design_matrix['TFL_DEATH_DT'])

	# set organ type to kidney only
	design_matrix = design_matrix[design_matrix.ORG_TY == "KI: Kidney"]
	design_matrix = design_matrix[design_matrix.REC_TX_ORG_TY == "KI: Kidney"]
		
	# recipient age is already includes only greater than 18
	# print(design_matrix['REC_AGE_AT_TX'].unique())
	# ['35-49' '50-64' '18-34' '65+']
	# recipient age over 18 - 216 months
	design_matrix = design_matrix[design_matrix['REC_AGE_IN_MONTHS_AT_TX'].astype(float) >= 216.0]
	print ("Recipient Age >=18: " + str(len(design_matrix)))

	# restrict donor age to greater than 9
	design_matrix = design_matrix[design_matrix.DON_AGE >= 9]
	print ("Donor Age >=9: " + str(len(design_matrix)))

	# Get death date.  Per Ann and Shannon (May 16, 2006) use pers_all_death_dt for all these purposes;
	# No longer have this, so use TFL_DEATH_DT;
	# Using additional follow-up (ESRD and SSDMF), so do not include tfl_lafudate, tfl_lafudateki;

	# print(design_matrix['TFL_DEATH_DT'].dtypes)

	design_matrix['TFL_DEATH_DT'] = pd.to_datetime(design_matrix['TFL_DEATH_DT'])
	design_matrix['TFL_GRAFT_DT'] = pd.to_datetime(design_matrix['TFL_GRAFT_DT'])
	design_matrix['PERS_RETX'] = pd.to_datetime(design_matrix['PERS_RETX'])
	design_matrix['TFL_ENDTXFU'] = pd.to_datetime(design_matrix['TFL_ENDTXFU'])
	design_matrix['REC_TX_DT'] = pd.to_datetime(design_matrix['REC_TX_DT'])

	# TODO - review censoring date - selected as end of year 2022 for pubsaf2306
	censor_dt = datetime.strptime("2022-12-31",'%Y-%m-%d')
	design_matrix['censor_dt'] = censor_dt
	design_matrix['died_dt'] = design_matrix['TFL_DEATH_DT']
	design_matrix.loc[ design_matrix['died_dt'] <= design_matrix['REC_TX_DT'], 'died_dt'] = np.nan # Death before the tx date is wrong;

	# if died_dt = < REC_TX_DT then died_dt = .; 
	design_matrix['endpt'] = design_matrix[['TFL_ENDTXFU','PERS_RETX','censor_dt','died_dt']].min(axis=1)
	design_matrix['patyrs'] = (design_matrix['endpt'] - design_matrix['REC_TX_DT']) / np.timedelta64(1, 'Y')

	# SAS code created died variable for recipient - yes or no variable
	# died=(died_dt=endpt) ;

	design_matrix['graftenddate'] = design_matrix[['TFL_GRAFT_DT','died_dt','TFL_ENDTXFU','PERS_RETX','censor_dt']].min(axis=1)
	design_matrix['graftyrs'] = (design_matrix['graftenddate'] - design_matrix['REC_TX_DT']) / np.timedelta64(1, 'Y')

	design_matrix['grf_fail'] = 0

	# only measure graft failure in first year
	# if (graftenddate=TFL_GRAFT_DT or graftenddate=died_dt or graftenddate=PERS_RETX) then grf_fail=1;  else grf_fail=0;
	design_matrix.loc[((design_matrix['graftenddate'] == design_matrix['TFL_GRAFT_DT']) | (design_matrix['graftenddate'] == design_matrix['died_dt']) | (design_matrix['graftenddate'] == design_matrix['PERS_RETX']) | (design_matrix['graftyrs'] <= 1)), 'grf_fail'] = 1


	# restrict transplant date
	design_matrix = design_matrix[design_matrix.REC_TX_DT.dt.year >= 2005] 
	print ("Transplant Date starting 2005: " + str(len(design_matrix)))

	# first transplants only
	design_matrix = design_matrix[design_matrix['REC_PREV_KI'] == 0] 
	print ("Transplant Date starting 2005: " + str(len(design_matrix)))
	design_matrix = design_matrix[design_matrix['REC_PREV_KP'] == 0]   # no changes

	# construct donor variables

	# %let ptkidonorvars = 		shared 		dcd		don_age donage_slope_ge18		dcadcodanox dcadcodcva dcadcodcnst dcadcodoth dcadcodhead		don_cmv_negative 		don_htn_0c
	#		ln_don_wgt_kg_0c ln_don_wgt_kg_0c_s55 don_wgt_kg_m		don_ecd		age_ecd		living_donor living_unrelated living_relate_missing living_related;

	# restrict donor race to White
	#design_matrix = design_matrix[design_matrix.DON_RACE == "8: White"]

	# restrict donor type to 'C' - already done
	# print(design_matrix['DON_TY'].unique())
	# ['C']

	# No relationship type because all deceased donors
	#  print(design_matrix['DON_RELATIONSHIP_TY'].unique())
	# [nan]

	# construct recipient variables
	# %let ptkirecvars = yearslice rec_age_at_tx rec_age_spline_35 rec_age_spline_50 rec_age_spline_65
	#                    diab_noted age_diab dm_can_age_spline_50
	# 				   can_dgn_htn_ndm can_dgn_pk_ndm can_dgn_gd_ndm
	#                    rec_prev_ki_tx rec_prev_ki_tx_dm
	#                    rbmi_0c rbmi_miss rbmi_gt_20 rbmi_DM rbmi_gt_20_dm
	#                    ln_c_hd_m ln_c_hd_0c ln_c_hd_m_ptx PKPRA_MS PKPRA_1080 PKPRA_GE80 PKPRA_lt10
	#                    hispanic CAN_RACE_BLACK CAN_RACE_asian CAN_RACE_WHITE;

	# restrict recipient race to White
	#sdesign_matrix = design_matrix[design_matrix.CAN_RACE == "8: White"] # dropped from 79238 to 36815

	# construct more HLA variables
	# %let ptkihlavars = mm0 mmDR0 mmDR1 mmA0 mmA1 mmB0 mmB1;

	# drop unneeded ID variables
	design_matrix = design_matrix.drop(['TX_ID'], axis=1)
	design_matrix = design_matrix.drop(['PX_ID'], axis=1)
	design_matrix = design_matrix.drop(['PERS_ID'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_ID'], axis=1)


	# drop unneeded HLA variables
	design_matrix = design_matrix.drop(['REC_A_MM_EQUIV_CUR'], axis=1)
	design_matrix = design_matrix.drop(['REC_B_MM_EQUIV_CUR'], axis=1)
	design_matrix = design_matrix.drop(['REC_DR_MM_EQUIV_CUR'], axis=1)
	design_matrix = design_matrix.drop(['HAPPAIR_RECIP'], axis=1)
	design_matrix = design_matrix.drop(['HAPPAIR_DONOR'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_HAP1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_HAP2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_HAP1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_HAP2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_A_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_A_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_A_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_A_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_B_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_B_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_B_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_B_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_C_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_C_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_C_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_C_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DRB1_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DRB1_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DRB1_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DRB1_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DRW1_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DRW1_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DRW1_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DRW1_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DQA1_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DQA1_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DQA1_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DQA1_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DQB1_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DQB1_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DQB1_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DQB1_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DPA1_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DPA1_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DPA1_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DPA1_2'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DPB1_1'], axis=1)
	design_matrix = design_matrix.drop(['RECIP_DPB1_2'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DPB1_1'], axis=1)
	design_matrix = design_matrix.drop(['DONOR_DPB1_2'], axis=1)

	# write SRTR design matrix
	design_matrix_filename = "SRTR_AA_MM_9loc_grffail_" + str(i) + ".txt"
	design_matrix.to_csv(design_matrix_filename,index=False,sep="\t")
	print("Logging info")
	
	#file.close()
	#exit() # TEMP - stop after one runmatch file
