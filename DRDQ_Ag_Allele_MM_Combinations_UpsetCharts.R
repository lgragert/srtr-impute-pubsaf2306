###############################################################################
#   SCRIPT NAME:    DRDQ_Ag_Allele_MM_Combinations_UpsetCharts.R
#   DESCRIPTION:    Enumerate SRTR Ag and Allele MM for DR, DQ, DQA1 Combinations
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
SRTR$DR_0Ag <- ifelse(SRTR$DR_AgMM == 0, 1, 0)
SRTR$DR_1Ag <- ifelse(SRTR$DR_AgMM == 1, 1, 0)
SRTR$DR_2Ag <- ifelse(SRTR$DR_AgMM == 2, 1, 0)
SRTR$DQ_0Ag <- ifelse(SRTR$DQ_AgMM == 0, 1, 0)
SRTR$DQ_1Ag <- ifelse(SRTR$DQ_AgMM == 1, 1, 0)
SRTR$DQ_2Ag <- ifelse(SRTR$DQ_AgMM == 2, 1, 0)
SRTR$DQA1_0Ag <- ifelse(SRTR$DQA1_AgMM == 0, 1, 0)
SRTR$DQA1_1Ag <- ifelse(SRTR$DQA1_AgMM == 1, 1, 0)
SRTR$DQA1_2Ag <- ifelse(SRTR$DQA1_AgMM == 2, 1, 0)

# Now for Allele MM
SRTR$DR_0Allele <- ifelse(SRTR$DR_AlleleMM == 0, 1, 0)
SRTR$DR_1Allele <- ifelse(SRTR$DR_AlleleMM == 1, 1, 0)
SRTR$DR_2Allele <- ifelse(SRTR$DR_AlleleMM == 2, 1, 0)
SRTR$DQ_0Allele <- ifelse(SRTR$DQ_AlleleMM == 0, 1, 0)
SRTR$DQ_1Allele <- ifelse(SRTR$DQ_AlleleMM == 1, 1, 0)
SRTR$DQ_2Allele <- ifelse(SRTR$DQ_AlleleMM == 2, 1, 0)
SRTR$DQA1_0Allele <- ifelse(SRTR$DQA1_AlleleMM == 0, 1, 0)
SRTR$DQA1_1Allele <- ifelse(SRTR$DQA1_AlleleMM == 1, 1, 0)
SRTR$DQA1_2Allele <- ifelse(SRTR$DQA1_AlleleMM == 2, 1, 0)

head(SRTR)

# 0-DRMM vs +1 MM
SRTR[SRTR == 2] <- 1

head(SRTR)

# save most plots independently with width=2000, (height varies, start with height=700)
# All onto one chart for all MM, 2 MM accounted for as 1 count
UpSetR::upset(SRTR, sets = c("DR_AgMM","DR_AlleleMM","DQA1_AgMM","DQA1_AlleleMM","DQ_AgMM","DQ_AlleleMM"), 
              sets.bar.color = c("magenta", "red", "orange","green","blue", "purple"),	
              point.size = 3.5, line.size = 2, intersection_matrix(geom=geom_point(shape='fill'), size=3.5,stroke=0.45),
              mainbar.y.label = "Instances of Ag-MM and Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: Allele_Ag_MM_DR-DQ-DQA1_Upset.png, width=2000, height=1500



# All cases for 0,1,2 into one plot
UpSetR::upset(SRTR, sets = c("DR_0Ag","DR_0Allele","DQA1_0Ag","DQA1_0Allele", "DQ_0Ag","DQ_0Allele", 
                             "DR_1Ag","DR_1Allele","DQA1_1Ag","DQA1_1Allele","DQ_1Ag","DQ_1Allele", 
                             "DR_2Ag","DR_2Allele","DQA1_2Ag","DQA1_2Allele","DQ_2Ag","DQ_2Allele"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink",
                                 "navy", "brown", "violet", "gray", 'gold', 'darkgreen', "maroon", "sienna", "skyblue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM and Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: 0-1-2_DR-DQ-DQA1_Ag_Allele_MM_Upset.PNG, width=5000, height=2000


# All cases for 0,1,2 but no DQA1
UpSetR::upset(SRTR, sets = c("DR_0Ag","DR_0Allele", "DR_1Ag","DR_1Allele","DR_2Ag","DR_2Allele",
                             "DQ_0Ag","DQ_0Allele", "DQ_1Ag","DQ_1Allele", "DQ_2Ag","DQ_2Allele"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink",
                                 "navy", "brown", "violet"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM and Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: DR-DQ_Allele_AgMM_Upset.PNG, width=2000, height=1500


# Only Ag Cases
UpSetR::upset(SRTR, sets = c("DR_0Ag","DR_1Ag","DR_2Ag","DQA1_0Ag", "DQA1_1Ag","DQA1_2Ag","DQ_0Ag","DQ_1Ag","DQ_2Ag"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM", sets.x.label = "Total Ag-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: 0-1-2_DR-DQA1-DQ_AgMM_Only_Upset.PNG


# No DQA1 AgMM
UpSetR::upset(SRTR, sets = c("DR_0Ag","DR_1Ag","DR_2Ag","DQ_0Ag","DQ_1Ag","DQ_2Ag"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM", sets.x.label = "Total Ag-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: No_DQA1_0-1-2_DR-DQ_AgMM_Upset.PNG


# Only Allele Cases
UpSetR::upset(SRTR, sets = c("DR_0Allele", "DR_1Allele","DR_2Allele","DQA1_0Allele","DQA1_1Allele","DQA1_2Allele","DQ_0Allele","DQ_1Allele", "DQ_2Allele"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Allele-MM", sets.x.label = "Total Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: 0-1-2_DR-DQA1-DQ_AlleleMM_Only_Upset.PNG


# No DQA1 Allele MM
UpSetR::upset(SRTR, sets = c("DR_0Allele", "DR_1Allele","DR_2Allele","DQ_0Allele","DQ_1Allele", "DQ_2Allele"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Allele-MM", sets.x.label = "Total  Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: No_DQA1_0-1-2_DR-DQ_AlleleMM_Upset.PNG



# Only 0-DR AgMM but with all other cases
UpSetR::upset(SRTR, sets = c("DQA1_0Ag","DQA1_1Ag","DQA1_2Ag","DQA1_0Allele", 
                             "DQA1_1Allele","DQA1_2Allele","DQ_0Ag","DQ_1Ag","DQ_2Ag","DQ_0Allele", 
                             "DQ_1Allele", "DQ_2Allele", "DR_0Allele","DR_1Allele","DR_2Allele","DR_0Ag"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink",
                                 "navy", "brown", "violet", 'gold', 'darkgreen', "maroon", "skyblue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM and Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: 0-DR_AgMM_all-cases_Upset.PNG



# Only 1-DR AgMM but with all other cases
UpSetR::upset(SRTR, sets = c("DQA1_0Ag","DQA1_1Ag","DQA1_2Ag","DQA1_0Allele", 
                             "DQA1_1Allele","DQA1_2Allele","DQ_0Ag","DQ_1Ag","DQ_2Ag","DQ_0Allele", 
                             "DQ_1Allele", "DQ_2Allele", "DR_0Allele","DR_1Allele","DR_2Allele","DR_1Ag"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink",
                                 "navy", "brown", "violet", 'gold', 'darkgreen', "maroon", "skyblue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM and Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: 1-DR_AgMM_all-cases_Upset.PNG



# Only 2-DR AgMM but with all other cases
UpSetR::upset(SRTR, sets = c("DQA1_0Ag","DQA1_1Ag","DQA1_2Ag","DQA1_0Allele", 
                             "DQA1_1Allele","DQA1_2Allele","DQ_0Ag","DQ_1Ag","DQ_2Ag","DQ_0Allele", 
                             "DQ_1Allele", "DQ_2Allele", "DR_0Allele","DR_1Allele","DR_2Allele","DR_2Ag"), 
              sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink",
                                 "navy", "brown", "violet", 'gold', 'darkgreen', "maroon", "skyblue"),	
              point.size = 3.5, line.size = 2,
              mainbar.y.label = "Instances of Ag-MM and Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
              text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# File name: 2-DR_AgMM_all-cases_Upset.PNG




# The ones that didn't make it, but helped me get to the ones above
# 0-DR AgMM with (0,1) DQB1, DQA1 AgMMs
# UpSetR::upset(SRTR, sets = c("DQA1_AgMM", "DQ_AgMM", "DR_0AgMM"), sets.bar.color = c("red","yellow","blue"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 0-DR Ag-MM with DQ Ag-MM", sets.x.label = "Total Ag-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# # save independently as: 0-DR_Ag-MM_DQ_Ag-MM_Upset.png, width=2000, height=700
# 
# 
# # 0-DR AgMM with (0,1) DQB1, DQA1 Allele MMs
# UpSetR::upset(SRTR, sets = c("DR_0AgMM","DQA1_AlleleMM", "DQ_AlleleMM"), sets.bar.color = c("red","yellow","blue"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 0-DR Ag-MM with DQ Allele-MM", sets.x.label = "Total 0 Ag-MM or Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# # save independently as: 0-DR_Ag-MM_DQ_Allele-MM_Upset.png width=2000, height=700
# 
# 
# # 0-DR Allele MM with (0,1) DQB1, DQA1 AgMMs
# UpSetR::upset(SRTR, sets = c("DR_0AlleleMM","DQA1_AgMM", "DQ_AgMM"), sets.bar.color = c("red","yellow","blue"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 0-DR Allele-MM with DQ Ag-MM", sets.x.label = "Total 0 Allee-MM or Ag-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# 
# # 0-DR Allele MM with (0,1) DQB1, DQA1 Allele MMs
# UpSetR::upset(SRTR, sets = c("DR_0AlleleMM","DQA1_AlleleMM","DQ_AlleleMM"), sets.bar.color = c("red","yellow","blue"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 0-DR Allele-MM with DQ Allele-MM", sets.x.label = "Total Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
#
#
# # Create plots for 0-DR Ag-MM agaisnt each 0,1,2 Ag-MM of DQ's
# UpSetR::upset(SRTR, sets = c("DR_0AgMM", "DR_1AgMM", "DR_2AgMM", "DQA1_0AgMM", "DQA1_1AgMM", "DQA1_2AgMM", "DQ_0AgMM", "DQ_1AgMM", "DQ_2AgMM"), 
#               sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of Ag-MM", sets.x.label = "Total Ag-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# 
# # 0-DR Ag-MM against each 0,1,2 Allele-MM of DQ's
# UpSetR::upset(SRTR, sets = c("DR_0AgMM", "DR_1AgMM", "DR_2AgMM", "DQ_0AlleleMM", "DQ_1AlleleMM", "DQ_2AlleleMM", "DQA1_0AlleleMM", "DQA1_1AlleleMM", "DQA1_2AlleleMM"), 
#               sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of DR Ag-MM or DQ's Allele-MM", sets.x.label = "Total Ag-MM or Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# # 0-DR Allele-MM against each 0,1,2 Ag-MM of DQ's
# UpSetR::upset(SRTR, sets = c("DR_0AlleleMM", "DR_1AlleleMM", "DR_2AlleleMM", "DQA1_0AgMM", "DQA1_1AgMM", "DQA1_2AgMM", "DQ_0AgMM", "DQ_1AgMM", "DQ_2AgMM"), 
#               sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of Allele-MM or Ag-MM", sets.x.label = "Total Allele-MM or Ag-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# 
# # 0-DR Allele-MM against each 0,1,2 Allele-MM of DQ's
# UpSetR::upset(SRTR, sets = c("DR_0AlleleMM", "DR_1AlleleMM", "DR_2AlleleMM", "DQ_0AlleleMM", "DQ_1AlleleMM", "DQ_2AlleleMM", "DQA1_0AlleleMM", "DQA1_1AlleleMM", "DQA1_2AlleleMM"), 
#               sets.bar.color = c("magenta", "red", "orange", "yellow","green","blue", "purple", "cyan", "pink"), point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of Allele-MM", sets.x.label = "Total Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# 
# # All onto one chart for 0 MM
# UpSetR::upset(SRTR, sets = c("DR_0AgMM","DR_0AlleleMM","DQA1_0AgMM","DQA1_0AlleleMM","DQ_0AgMM","DQ_0AlleleMM"), 
#               sets.bar.color = c("magenta", "red", "orange","green","blue", "purple"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 0 Ag-MM and 0 Allele-MM", sets.x.label = "Total 0 Ag-MM or 0 Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# 
# # All onto one chart for 1 MM
# UpSetR::upset(SRTR, sets = c("DR_1AgMM","DR_1AlleleMM","DQA1_1AgMM","DQA1_1AlleleMM","DQ_1AgMM","DQ_1AlleleMM"), 
#               sets.bar.color = c("magenta", "red", "orange","green","blue", "purple"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 1 Ag-MM and 1 Allele-MM", sets.x.label = "Total 1 Ag-MM or 1 Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
# 
# 
# # All onto one chart for 2 MM
# UpSetR::upset(SRTR, sets = c("DR_2AgMM","DR_2AlleleMM","DQA1_2AgMM","DQA1_2AlleleMM","DQ_2AgMM","DQ_2AlleleMM"), 
#               sets.bar.color = c("magenta", "red", "orange","green","blue", "purple"),	
#               point.size = 3.5, line.size = 2,
#               mainbar.y.label = "Instances of 2 Ag-MM and 2 Allele-MM", sets.x.label = "Total 2 Ag-MM or 2 Allele-MM",
#               text.scale = c(1.3, 1.3, 1, 1, 2, 1.60), empty.intersections ="on", order.by = "freq",keep.order = TRUE)
