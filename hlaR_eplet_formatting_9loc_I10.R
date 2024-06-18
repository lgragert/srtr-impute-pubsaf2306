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
ARRAY_ID<-paste0(split <- command_args[1])
print(paste0("SLURMARRAYID: ",ARRAY_ID))

filename <- paste0("./SRTR_AA_MM_9loc_matrix_", ARRAY_ID, ".txt")
print(paste0("filename:",filename ))


#NEED to set values to specific row counts!
# myeach = number of rows you want in a table, must be even to include both D|R in analysis
# my rep and myseq = number of data frames formed and must be same value
#if not set correctly will give:   data length is not a multiple of split variable
#however it qill still include all rows one list may just be shorter
myeach = paste0(24756)
myrep = paste0(1:14)
myseq = paste0(1:14)

#Set up Matrix to correct format for analysis by functions
SRTR <- read.delim(filename, sep = '\t',header = T)

print(colnames(SRTR))

#Select Imputed Genotype of D|R and format for hlaR functions
#Columns selected in order for hlaR classII renaming
recip <- SRTR[, c("PX_ID","RECIP_DRB345_1","RECIP_DRB345_2","RECIP_DRB1_1","RECIP_DRB1_2",
                          "RECIP_DQA1_1","RECIP_DQA1_2","RECIP_DQB1_1","RECIP_DQB1_2",
                          "RECIP_DPA1_1","RECIP_DPA1_2","RECIP_DPB1_1","RECIP_DPB1_2")]

DONOR <- SRTR[, c("PX_ID","DONOR_DRB345_1","DONOR_DRB345_2","DONOR_DRB1_1","DONOR_DRB1_2",
                          "DONOR_DQA1_1","DONOR_DQA1_2","DONOR_DQB1_1","DONOR_DQB1_2",
                          "DONOR_DPA1_1","DONOR_DPA1_2","DONOR_DPB1_1","DONOR_DPB1_2")]

#Column names for hlaR classII processing
names <- c("pair_id","DRw1","DRw2","DRB1","DRB2",
           "DQA1","DQA2","DQB1","DQB2","DPA1","DPA2","DPB1","DPB2")
           

colnames(recip)<-names
colnames(DONOR)<-names


DONOR$subject_type <- "donor"

recip$subject_type <- "recipient"


SRTR <- rbind(DONOR,recip)

#column order for hlaR processing
col.order <- c("pair_id","subject_type","DRB1","DRB2","DRw1","DRw2","DQB1","DQB2",
               "DQA1","DQA2","DPB1","DPB2","DPA1","DPA2")

SRTR <- SRTR[ , col.order]


SRTR <- SRTR %>% arrange(pair_id)

classII <-SRTR %>%  dplyr::select("pair_id","subject_type",starts_with("D"))

classII <- classII %>% arrange(pair_id)
#write.table(classII,file="./hlaRepletFormatCII.csv", sep = ",", row.names = T)


classII <- classII %>% dplyr::na_if("DRBX*NNNN")
classII[is.na(classII)] <- "" 


#classI <-SRTR %>%  dplyr::select("pair_id","subject_type",starts_with("a"),
                                 #starts_with("b"),starts_with("c"))

#classI <- classI %>% arrange(pair_id)
#write.table(classI,file="./hlaRepletFormatCI.csv", sep = ",", row.names = T)

ls <- setNames(split(classII, rep(myrep, each = myeach)), paste("hlaRClassIIFormat", seq(myseq), sep = "_"))


#how to extract and save sublists with specific names as dataframes from a list
for (nameframe in names(ls)){
  print(nameframe)
  formatted<-as.data.frame(ls[[nameframe]])
  write.table(formatted,file=paste0("./Imp_",ARRAY_ID,"_",nameframe,".csv"), sep = ",", row.names = T)
}

