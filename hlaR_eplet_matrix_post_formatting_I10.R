###############################################################################
#   SCRIPT NAME:    hlaR_eplet_matrix_post_formating_I10.R
#   DESCRIPTION:    combines generated eplet matrix for all 10 Imputations
#   OUTPUT:         EpletMM_Matrix_%a.csv
#   DATE:           February 23, 2023
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
if (!require("cli"))
  installed.packages("cli")
if (!require("devtools"))
  installed.packages("devtools")
if (!require("tidyselect"))
  install.packages("tidyselect")
if (!require("readr"))
  install.packages("readr")
if (!require("reshape2"))
  install.packages("reshape2")

#Set path
#path<- paste0(getwd(),"/")


#Step 1 fectch all necessary eplet MM matrix files and filter for
#hlaR ouput only with PXIDs in FIBERS Analysis

#fetch hlaR files

#for specific Eplets
# Fetch command line argument for the imputation ID
command_args <- commandArgs(trailingOnly = TRUE)
IMP_ID<-paste0(split <- command_args[1])

# Define the file name pattern for the hlaR eplet-only matrix files
filename <- paste0("Imp_",IMP_ID,"_hlaR_eplet_only_matrix")
#filename <- paste0("Imp_1_hlaR_eplet_only_matrix")
print(filename)
hlaRDDCII <- dir(path="./hlaR_March23", pattern=paste0(filename), full.names=TRUE, recursive = TRUE)
print(hlaRDDCII)

# Read and combine all matching eplet matrix files into a single dataframe
eplets <- do.call(`rbind`,lapply(hlaRDDCII, read.delim, sep=",",header = TRUE, fill = TRUE))
#eplets <- read.table(filename, sep = ',')

eplets[is.na(eplets)] <- 0 # Replace NA values with 0

# Read hlaR global matrix file and select relevant columns
#join to DR DQ typing from hlaR matrix
hlaR <- paste0("./hlaR_9loc_global_matrix_",IMP_ID,".txt")
#hlaR <- paste0("./hlaR_replicates/hlaR_9loc_matrix_1.txt")
hlaR <- read.delim(hlaR, sep = '\t',header=TRUE)

hlaR <- hlaR %>% dplyr::select("PX_ID",starts_with("RECIP_DR"),starts_with("RECIP_DQ"),
                       starts_with("DONOR_DR"),starts_with("DONOR_DQ"))

# Joine eplets data with the hlaR data
hlaR <- left_join(eplets,hlaR)


# join to generated 
# Define the filename pattern
filename <- paste0("IMP_",IMP_ID,"_hlaR_Weibe_regression_input.csv")
#filename <- paste0("IMP_1_hlaR_Weibe_regression_input_")
hlaRDDCII <- dir(path="./hlaR_March23", pattern=paste0(filename), full.names=TRUE, recursive = TRUE)

# Read and combine all matching Weibe input files into a single dataframe
SRTR1 <- do.call(`rbind`,lapply(hlaRDDCII, read.delim,sep=",",header = TRUE, fill = TRUE))

# Join the combined hlaR data with the SRTR1 data
SRTR1 <- left_join(SRTR1,hlaR)
#SRTR1 <- left_join(SRTR1,eplets)

#print(nrow(SRTR1))

#summary(SRTR1)

# Output the combined data to a text file
E_filename <- paste0("./hlaR_March23/EpletMM_Matrix_",IMP_ID, ".txt") 
write.table(SRTR1,file=paste0(E_filename), sep = "\t", row.names = F)

#Runmatch format for Kieth

#debugging locally remove
# Read a reference SRTR file to obtain a list of PX_IDs
SRTR <- read.table("./SRTR_AA_MM_matrix_grffail_replicates_2022-04-03/SRTR_AA_MM_matrix_grffail_1.txt",sep = '\t',header = T)
id_list <-as.list(SRTR$PX_ID)

#debugging locally remove
#colnames(SRTR1)[which(names(SRTR1) == "risk.y")] <- "hlaR_risk"
# Filter eplets data for IDs present in the SRTR reference file and reshape data to long format
eplets<-filter(eplets, PX_ID %in% id_list)
eplets <- pivot_longer(eplets, cols = c(2:275), names_to = "Eplets", values_to = "Mismatches")

# Add columns for mismatch counts
eplets$mm_count_0 <- ''
eplets$mm_count_1 <- ''
eplets$mm_count_2 <- ''

for (i in 1:nrow(eplets)){
if(eplets$Mismatches[i]==0){
  eplets$mm_count_0[i]<- 1
  eplets$mm_count_1[i]<- 0
  eplets$mm_count_2[i]<- 0
}
  if(eplets$Mismatches[i]==1){
    eplets$mm_count_0[i] <- 0
    eplets$mm_count_1[i] <- 1
    eplets$mm_count_2[i] <- 0
    
  }
  if(eplets$Mismatches[i]==2){
    eplets$mm_count_0[i] <- 0
    eplets$mm_count_1[i] <- 0
    eplets$mm_count_2[i] <- 1
    
  }
}

print(summary(eplets))
runmatch <- eplets[,-3]
print(summary(runmatch))

# Output the runmatch data to a text file
E_filename <- paste0("./hlaR_March23/EpletMM_Matrix_RunMatch_",IMP_ID, ".txt") 
write.table(runmatch,file=paste0(E_filename), sep = "\t", row.names = F)
