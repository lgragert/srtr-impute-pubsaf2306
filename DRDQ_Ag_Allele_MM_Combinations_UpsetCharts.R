###############################################################################
#   SCRIPT NAME:    DRDQ_Ag_Allele_MM_Combinations_UpsetCharts.R
#   DESCRIPTION:    Enumerate SRTR 0-DR Ag and Allele MM with DQ Ag and Allele MM Combinations
#   OUTPUT:         *_Upset.png
#   DATE:           January 17, 2024
#   AUTHOR:         Grace Lord (gwager@tulane.edu ; GitHub: gwager), modified by Alyssa Paynter
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
if (!require("pastecs"))
  installed.packages("pastecs")
if (!require("tidyselect"))
  install.packages("tidyselect")
if (!require("ggforce"))
  install.packages("ggforce")
library(ggforce) 
if (!require(devtools)) 
  install.packages("devtools")
if (!require("ComplexUpset"))
  install.packages("ComplexUpset")
if (!require("UpSetR"))
  install.packages("UpSetR")

#Utilization of ComplexUpset 
#https://stackoverflow.com/questions/67094573/upsetr-sets-bar-interacting-different-color-and-its-sets-intersections
#https://github.com/krassowski/complex-upset/
#xhttps://krassowski.github.io/complex-upset/articles/Examples_R.html
#https://krassowski.github.io/complex-upset/reference/upset.html
#https://github.com/krassowski/complex-upset/issues


filename <- paste0("./srtr_ag_allele_mm_1.csv")
SRTR <- read.table(filename,sep = ',',header = T,quote='', comment='')

SRTR <- SRTR %>% dplyr:: select('PX_ID', 'REC_DR_MM_EQUIV_CUR', 'REC_DQ_MM_EQUIV_CUR', 'REC_DQA1_MM_EQUIV_CUR', 'REC_DR_ALLELE_MM', 'REC_DQ_ALLELE_MM', 'REC_DQA1_ALLELE_MM')

# Change the names of columns to make them easier to work with
names(SRTR)[names(SRTR) == 'REC_DR_MM_EQUIV_CUR'] <- 'DR_AgMM'
names(SRTR)[names(SRTR) == 'REC_DQ_MM_EQUIV_CUR'] <- 'DQ_AgMM'
names(SRTR)[names(SRTR) == 'REC_DQA1_MM_EQUIV_CUR'] <- 'DQA1_AgMM'
names(SRTR)[names(SRTR) == 'REC_DR_ALLELE_MM'] <- 'DR_AlleleMM'
names(SRTR)[names(SRTR) == 'REC_DQ_ALLELE_MM'] <- 'DQ_AlleleMM'
names(SRTR)[names(SRTR) == 'REC_DQA1_ALLELE_MM'] <- 'DQA1_AlleleMM'


# Make table data frame and make it all numeric
SRTR <-as.data.frame(lapply(SRTR,as.numeric))

# Create separate columns for each type of mismatch we have
SRTR$DR_0AgMM <- ifelse(SRTR$DR_AgMM == 0, 1, 0)
SRTR$DR_1AgMM <- ifelse(SRTR$DR_AgMM == 1, 1, 0)
SRTR$DR_2AgMM <- ifelse(SRTR$DR_AgMM == 2, 1, 0)
SRTR$DQ_0AgMM <- ifelse(SRTR$DQ_AgMM == 0, 1, 0)
SRTR$DQ_1AgMM <- ifelse(SRTR$DQ_AgMM == 1, 1, 0)
SRTR$DQ_2AgMM <- ifelse(SRTR$DQ_AgMM == 2, 1, 0)
SRTR$DQA1_0AgMM <- ifelse(SRTR$DQA1_AgMM == 0, 1, 0)
SRTR$DQA1_1AgMM <- ifelse(SRTR$DQA1_AgMM == 1, 1, 0)
SRTR$DQA1_2AgMM <- ifelse(SRTR$DQA1_AgMM == 2, 1, 0)

# Now for Allele MM
SRTR$DR_0AlleleMM <- ifelse(SRTR$DR_AlleleMM == 0, 1, 0)
SRTR$DR_1AlleleMM <- ifelse(SRTR$DR_AlleleMM == 1, 1, 0)
SRTR$DR_2AlleleMM <- ifelse(SRTR$DR_AlleleMM == 2, 1, 0)
SRTR$DQ_0AlleleMM <- ifelse(SRTR$DQ_AlleleMM == 0, 1, 0)
SRTR$DQ_1AlleleMM <- ifelse(SRTR$DQ_AlleleMM == 1, 1, 0)
SRTR$DQ_2AlleleMM <- ifelse(SRTR$DQ_AlleleMM == 2, 1, 0)
SRTR$DQA1_0AlleleMM <- ifelse(SRTR$DQA1_AlleleMM == 0, 1, 0)
SRTR$DQA1_1AlleleMM <- ifelse(SRTR$DQA1_AlleleMM == 1, 1, 0)
SRTR$DQA1_2AlleleMM <- ifelse(SRTR$DQA1_AlleleMM == 2, 1, 0)

head(SRTR)

# 0-DRMM vs +1 MM
SRTR[SRTR == 2] <- 1

head(SRTR)

# Start with 0-DR MM vs +1 DQ MMs
# 0-DR AgMM with (0,1) DQB1, DQA1 AgMMs
UpSetR::upset(SRTR, sets = c("DQA1_AgMM", "DQ_AgMM", "DR_0AgMM"), sets.bar.color = c("red","yellow","blue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of 0-DR Ag-MM with DQ Ag-MM", sets.x.label = "Total Ag-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# save independently as: 0-DR_Ag-MM_DQ_Ag-MM_Upset.png, width=2000, height=700


# 0-DR AgMM with (0,1) DQB1, DQA1 Allele MMs
UpSetR::upset(SRTR, sets = c("DR_0AgMM","DQA1_AlleleMM", "DQ_AlleleMM"), sets.bar.color = c("red","yellow","blue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of 0-DR Ag-MM with DQ Allele-MM", sets.x.label = "Total 0 Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# save independently as: 0-DR_Ag-MM_DQ_Allele-MM_Upset.png width=2000, height=700


# 0-DR Allele MM with (0,1) DQB1, DQA1 AgMMs
UpSetR::upset(SRTR, sets = c("DR_0AlleleMM","DQA1_AgMM", "DQ_AgMM"), sets.bar.color = c("red","yellow","blue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of 0-DR Allele-MM with DQ Ag-MM", sets.x.label = "Total 0 Allee-MM or Ag-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)


# 0-DR Allele MM with (0,1) DQB1, DQA1 Allele MMs
UpSetR::upset(SRTR, sets = c("DR_0AlleleMM","DQA1_AlleleMM","DQ_AlleleMM"), sets.bar.color = c("red","yellow","blue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of 0-DR Allele-MM with DQ Allele-MM", sets.x.label = "Total Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)





# Create plots for 0-DR Ag-MM agaisnt each 0,1,2 Ag-MM of DQ's
UpSetR::upset(SRTR, sets = c("DR_0AgMM", "DR_1AgMM", "DR_2AgMM", "DQA1_0AgMM", "DQA1_1AgMM", "DQA1_2AgMM", "DQ_0AgMM", "DQ_1AgMM", "DQ_2AgMM"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM", sets.x.label = "Total Ag-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)


# 0-DR Ag-MM against each 0,1,2 Allele-MM of DQ's
UpSetR::upset(SRTR, sets = c("DR_0AgMM", "DR_1AgMM", "DR_2AgMM", "DQ_0AlleleMM", "DQ_1AlleleMM", "DQ_2AlleleMM", "DQA1_0AlleleMM", "DQA1_1AlleleMM", "DQA1_2AlleleMM"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of DR Ag-MM or DQ's Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)

# 0-DR Allele-MM against each 0,1,2 Ag-MM of DQ's
UpSetR::upset(SRTR, sets = c("DR_0AlleleMM", "DR_1AlleleMM", "DR_2AlleleMM", "DQA1_0AgMM", "DQA1_1AgMM", "DQA1_2AgMM", "DQ_0AgMM", "DQ_1AgMM", "DQ_2AgMM"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Allele-MM or Ag-MM", sets.x.label = "Total Allele-MM or Ag-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)


# 0-DR Allele-MM against each 0,1,2 Allele-MM of DQ's
UpSetR::upset(SRTR, sets = c("DR_0AlleleMM", "DR_1AlleleMM", "DR_2AlleleMM", "DQ_0AlleleMM", "DQ_1AlleleMM", "DQ_2AlleleMM", "DQA1_0AlleleMM", "DQA1_1AlleleMM", "DQA1_2AlleleMM"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Allele-MM", sets.x.label = "Total Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)




