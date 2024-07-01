###############################################################################
#   SCRIPT NAME:    hlaR_DRDQ_single_molecule_post.R
#   DESCRIPTION:    Extract hlaR Risk column as well as assign 
#                   Wiebe et al Risk Categories
#   OUTPUT:         hlaR_Weibe_regression_input.csv
#                   hlaR_Wiebe_FIBERS_input.csv
#   DATE:           January 31, 2023
#   AUTHOR:         Grace Lord (gwager@tulane.edu ; GitHub: gwager)
#   PI:             Loren Gragert, Ph.D.
#   ORGANIZATION:   Tulane University School of Medicine
#   NOTES:
###############################################################################

if(!require("dtplyr"))
  install.packages("dtplyr")
if(!require("dplyr"))
  install.packages("dplyr")
if (!require("tidyverse"))
  install.packages("tidyverse")
if (!require("dplyr"))
  install.packages("dplyr")
if (!require("plyr"))
  install.packages("plyr")
if (!require("tidyr"))
  install.packages("tidyr")
if (!require("stringr"))
  install.packages("stringr")
if (!require("scales"))
  install.packages("scales")
if (!require("data.table"))
  installed.packages("data.table")
if (!require("rlang"))
  installed.packages("rlang")
if (!require("tidyselect"))
  install.packages("tidyselect")
if (!require("Hmisc"))
  install.packages("Hmisc")
if (!require("mlbench"))
  install.packages("mlbench")
if (!require("caret"))
  install.packages("caret")
if (!require("leaps"))
  install.packages("leaps")
if (!require("tools"))
  install.packages("tools")
if (!require("purrr"))
  install.packages("purrr")
if (!require("fs"))
  install.packages("fs")
if (!require("matrixStats"))
  install.packages("matrixStats")
if (!require("modeest"))
  install.packages("modeest")
if (!require("reshape2"))
  install.packages("reshape2")
if (!require("questionr"))
  install.packages("questionr")
if (!require("nlme"))
  install.packages("nlme")


#Set path
path<- paste0(getwd(),"/")

#fetch hlaR files
#fetch hlaR files
command_args <- commandArgs(trailingOnly = TRUE)
IMP_ID<-paste0(split <- command_args[1])
Format_ID<-paste0(split <- command_args[2])

filename <- paste0("./hlaR_February_22_2023/Imp_",IMP_ID , "_hlaRepletsDRDQRiskA_", Format_ID, ".csv")
print(paste0("filename:",filename ))

molecule <- read.table(filename, sep = ',',header = TRUE)



#extract hlaR risk assignment and generate boolean columns for each
#   risk = ifelse(between(DQ, 15, 31), "high",
#   ifelse((DR >= 7 & DQ <= 14) | (DR < 7 & between(DQ, 9, 15)), "interm",
#   ifelse(DR < 7 & DQ < 9, "low", "out of bound"))))

hlaRisk <- molecule %>% dplyr::select("pair_id","risk")

hlaRisk$hlaR_Low <- 0
hlaRisk$hlaR_Medium <- 0
hlaRisk$hlaR_High <- 0

for (i in 1:nrow(hlaRisk)){
  if(hlaRisk$risk[i]== "low"){
    hlaRisk$hlaR_Low[i] <- 1
    
  }
  else if(hlaRisk$risk[i]== "interm"){
    hlaRisk$hlaR_Medium[i] <- 1
    
  }
  else if(hlaRisk$risk[i]== "high"){
    hlaRisk$hlaR_High[i] <- 1
    
  }else{
    name <- hlaRisk$pair_id[i]
    print(name)
  }
  
}


#low/med/high risk looping per wiebe et al. AJT 2019
# Figure 3.C values found to be correlated with DSAs and survival
molecule$Wiebe_Low <- 0
molecule$Wiebe_Medium <- 0
molecule$Wiebe_High <- 0

for (i in 1:nrow(molecule)){
  if(molecule$DR[i]==0){
    if(molecule$DQ[i]==0){
      molecule$Wiebe_Low[i] <- 1 
    }else if(molecule$DQ[i]>=1 & molecule$DQ[i]<=8){
      molecule$Wiebe_Medium[i] <- 1
    }else if (molecule$DQ[i]>=9) {
      molecule$Wiebe_High[i] <- 1
    }
  }
  if(molecule$DR[i]>=1 & molecule$DR[i]<=6){
    if(molecule$DQ[i]<=8){
      molecule$Wiebe_Medium[i] <- 1
    }else if(molecule$DQ[i]>=9) {
      molecule$Wiebe_High[i] <- 1
    }
  }
  if(molecule$DR[i]>=7){
    molecule$Wiebe_High[i] <- 1  }
}


#generate risk_category column for factor levels
molecule$wiebe_risk <- 0

for (i in 1:nrow(molecule)){
  if(molecule$Wiebe_Low[i]==1){
      molecule$wiebe_risk[i] <- "low"
    }else if(molecule$Wiebe_Medium[i]==1){
      molecule$wiebe_risk[i] <- "interm"
    }else if(molecule$Wiebe_High[i]==1){
      molecule$wiebe_risk[i] <- "high"
    }
}


WhlaR <- left_join(molecule,hlaRisk,by= "pair_id")
WhlaR <- unique(WhlaR)
colnames(WhlaR)[which(names(WhlaR) == "pair_id")] <- "PX_ID"
colnames(WhlaR)[which(names(WhlaR) == "risk")] <- "hlaR_risk"
summary(WhlaR)

E_filename <- paste0("./hlaR_February_22_2023/hlaR_Weibe_regression_input_",IMP_ID, ".csv") 
write.table(WhlaR,file=paste0(E_filename), sep = ",", row.names = F)

#WhlaR = subset(WhlaR, select = -c(DR,DQ) )

#SRTR <- read.table("./FIBERS_RISK_AgMM.csv", sep = ',',header = T)
#summary(SRTR)


#using this file as the pipeline was developed out of the order it is documented to flow
SRTR <- read.table("./KM_hlaR_Wiebe_FIBERS_input_1.txt", sep = ',',header = T)


WhlaR <- left_join(SRTR, WhlaR)
summary(WhlaR)
E_filename <- paste0("./hlaR_February_22_2023/KM_hlaR_Wiebe_FIBERS_",IMP_ID, ".txt") 
write.table(SRTR1,file=paste0(E_filename), sep = "\t", row.names = F)
