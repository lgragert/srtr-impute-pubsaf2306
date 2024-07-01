###############################################################################
#   SCRIPT NAME:    hlaR_eplet_IDA_I10.R
#   DESCRIPTION:    Process formatted Matrix for hlaR analysis on cypress
#   SLURM:          hlaR_eplet_IDA_I10_RUN.sh
#   OUTPUT:         ./IMP_*_hlaRepletsSumsclassIIA_*.csv
#   OUTPUT:         ./IMP_*_hlaRepletsListclassIIA_*.csv
#   OUTPUT:         ./IMP_*hlaRepletsDRDQRiskA_*.csv
#   DATE:           January 27, 2023
#   AUTHOR:         Grace Lord (gwager@tulane.edu ; GitHub: gwager)
#   PI:             Loren Gragert, Ph.D.
#   ORGANIZATION:   Tulane University School of Medicine
#   NOTES:
###############################################################################
#With SlURM script and setup as is 02.23.23 
#Can run on Cyp to alter all 10 imputations
#across 14 files within 



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
if (!require("cli"))
  installed.packages("cli")
if (!require("devtools"))
  installed.packages("devtools")
if (!require("tidyselect"))
  install.packages("tidyselect")
if (!require("readr"))
  install.packages("readr")
if (!require("hlaR"))
  install.packages("hlaR")

command_args <- commandArgs(trailingOnly = TRUE)
IMP_ID<-paste0(split <- command_args[1])
Format_ID<-paste0(split <- command_args[2])

filename <- paste0("./Imp_",IMP_ID,"_hlaRClassIIFormat_", Format_ID, ".csv")
print(paste0("filename:",filename ))

classII <- read.table(filename,sep = ',')



epiclassIIsingle <- data.frame(pair_id=numeric(), subject=numeric(), 
                              hla=character(),mm_eplets=character(),
                              mm_count=numeric())
epiclassIIoverall <- data.frame(pair_id=numeric(), subject=numeric(), 
                               mm_count_tt=numeric(), mm_count_unique=numeric())
epiclassIIrisk <- data.frame(pair_id=numeric(), DQ=numeric(), 
                             DR=numeric(), risk=character())


pair_id_list <- as.list(classII$pair_id)
pair_id_list <- unique (pair_id_list)

count <- (nrow(classII)/2)
count <-as.numeric(count)
print(paste0("N D|R pairs:", count))

  for (ID in pair_id_list){
    flag_ID <- paste0(ID)
    #print(ID)
    #made it an integer literal need to pull the L from end
    ID <-as.numeric(ID)
    flag_ID <- paste0(ID)
    
    
    #extract and process one set of IDs at a time
    PID <- filter(classII, pair_id == ID)
    PID<- as.data.frame(PID)
    PID <-CalEpletMHCII(PID, ver= 3)
    epirunSdf<-as.data.frame(PID[["single_detail"]])
    epirunOdf<-as.data.frame(PID[["overall_count"]])
    epirunRdf<-as.data.frame(PID[["dqdr_risk"]])
    epiclassIIoverall <- rbind(epiclassIIoverall,epirunOdf)
    epiclassIIsingle <- rbind(epiclassIIsingle,epirunSdf)
    epiclassIIrisk <- rbind(epiclassIIrisk,epirunRdf)
}

  
OA_filename <- paste0("./hlaR_February_22_2023/Imp_", IMP_ID,"_hlaRepletsSumsclassIIA_",Format_ID,".csv")
S_filename <- paste0("./hlaR_February_22_2023/Imp_", IMP_ID,"_hlaRepletsListclassIIA_",Format_ID ,".csv")
R_filename <- paste0("./hlaR_February_22_2023/Imp_", IMP_ID,"_hlaRepletsDRDQRiskA_",Format_ID,".csv")


write.table(epiclassIIsingle,file=paste0(S_filename), sep = ",", row.names = T)
write.table(epiclassIIoverall,file=paste0(OA_filename), sep = ",", row.names = T)
write.table(epiclassIIrisk,file=paste0(R_filename), sep = ",", row.names = T)
