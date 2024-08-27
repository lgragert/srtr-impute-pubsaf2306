###############################################################################
#   SCRIPT NAME:    Ethnic_Disparity.R
#   DESCRIPTION:    Plot mm count distribution for different genetic reference pops
#   OUTPUT:         DRDQ pdfs of freqs in wiebe and hlaR mm categories and FIBERS
#   DATE:           February 03, 2023
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
if (!require("pastecs"))
  installed.packages("pastecs")
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
if (!require("survival"))
  install.packages("survival")
if (!require("ranger"))
  install.packages("ranger")
library(ggfortify)

#Set path
path<- paste0(getwd(),"/")


#set risk categories to analyze
ethnicity <- c("16: Black or African American","64: Asian","8: White",
               "2000: Hispanic/Latino","128: Native Hawaiian or Other Pacific Islander",
               "Multi-Racial","32: American Indian or Alaska Native")

#fetch SRTR gf file
molecule <- read.table("./KM_hlaR_Wiebe_FIBERS_input_1.txt", sep = '\t',header = T)
risk_category_list <- molecule %>% dplyr::select("hlaR_risk", "wiebe_risk")

#https://stackoverfLow.com/questions/3872070/how-to-force-r-to-use-a-specified-factor-level-as-reference-in-a-regression
levels(molecule$risk_category)
levels(molecule$risk_category)[1]
molecule$wiebe_risk <- factor(molecule$wiebe_risk, levels = c("low","interm","high"))
molecule$hlaR_risk <- factor(molecule$hlaR_risk, levels = c("low","interm","high"))
molecule$I1_risk <- factor(molecule$I1_risk, levels = c("low","high"))
molecule$I9_risk <- factor(molecule$I9_risk, levels = c("low","high"))
molecule$g1g2_risk <- factor(molecule$g1g2_risk, levels = c("low","interm","high"))
molecule$REC_DR_MM_EQUIV_CUR <- factor(molecule$REC_DR_MM_EQUIV_CUR, levels = c(0,1,2))
molecule$REC_DQ_MM_EQUIV_CUR <- factor(molecule$REC_DQ_MM_EQUIV_CUR, levels = c(0,1,2))


# Bar plot of distribution of eplet mismatch counts for DR and for DQ.
Drc <- molecule %>% dplyr::select("PX_ID","DR", "CAN_RACE")
Drc$loci <- "DR"
names(Drc)[2] <- "mm_cnt"

DQc <- molecule %>% dplyr::select("PX_ID","DQ", "CAN_RACE")
DQc$loci <- "DQ"
names(DQc)[2] <- "mm_cnt"

 
DrDq <- rbind(Drc, DQc)


DrDq[DrDq == "16: Black or African American"]<- "AFA"
DrDq[DrDq == "64: Asian"]<- "ASN"
DrDq[DrDq == "8: White"]<- "CAU"
DrDq[DrDq == "2000: Hispanic/Latino"]<- "HIS"
DrDq[DrDq == "128: Native Hawaiian or Other Pacific Islander"]<- "HPI"
DrDq[DrDq == "Multi-Racial"]<- "MLT"
DrDq[DrDq == "32: American Indian or Alaska Native"]<- "NAM"

#barbell plot
DrDq %>% 
  ggplot(aes(x=CAN_RACE,
             y=mm_cnt,
             color=loci))+
  geom_boxplot(lwd=1)+
  #geom_jitter(width=0.1,alpha=0.2)+
  theme(legend.position = "right")+
  theme_bw()+
  xlab("Transplant Recipient Ethnicity")+
  ylab("Single Molecule Eplet Mismatch Count")+
  labs(title="DR|DQ Single Molecule Eplet Mismatch Counts 
       by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="DRDQ_MM_Counts_BP_up.pdf", width = 10, height = 10, units="in", path = "./")




#% of TX pairs in Risk Cat by ethnicity plot: Disparity interest
pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","DON_RACE","DR","DQ","wiebe_risk")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, wiebe_risk == "low")
  ln <- nrow(Low)
  pln <- ln/N*100
  med <- filter(pops, wiebe_risk == "interm")
  mn <- nrow(med)
  pmn <- mn/N*100
  High <- filter(pops, wiebe_risk == "high")
  hn <- nrow(High)
  phn <- hn/N*100

  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }


  Percentage <- c(pln,pmn,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("Low", "Medium", "High")
  pop_row_df$wiebe_risk <- risks
  pop_table<- rbind(pop_table, pop_row_df)

  Counts <- c(ln,mn,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("Low", "Medium", "High", "Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_wiebe_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_wiebe_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_wiebe_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)

#pop_table$wiebe_risk <- factor(pop_table$wiebe_risk, levels = c("High", "Medium", "Low"))

ggplot(data = pop_table, aes( x = wiebe_risk, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("Wiebe et Al. DR|DQ Single Molecule Eplet Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="Wiebe_DRDQ_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = wiebe_risk)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("RWiebe et Al. DR|DQ Single Molecule Eplet Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="Wiebe_DRDQ_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")



WN = 93+73+498
wpln= 93/WN*100
wpmn =73/WN*100
wphn = 498/WN*100

Percentage <- c(wpln,wpmn,wphn)
wiebe <- as.data.frame(Percentage)
wiebe$Ethnicity <- "Weibe et al. Cau"
risks <- c("Low", "Medium", "High")
wiebe$wiebe_risk <- risks
wiebe<- rbind(pop_table, wiebe)

#ensure order of risk cats
wiebe$wiebe_risk <- factor(wiebe$wiebe_risk, levels = c("High", "Medium", "Low"))


ggplot(data = wiebe, aes( x = wiebe_risk, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  ggtitle("Wiebe et Al. DR|DQ Single Molecule Eplet Risk Categories
              by Recipient Ethnic Population")+
  theme_bw()+
  xlab("risk_category")+
  ylab("Percent of Transplant Recipients")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="Wiebe_compare_all_clustered.pdf", width = 10, height = 10, units="in", path = "./")


L <- c("CAU","Weibe et al. Cau")

wiebe <- subset(wiebe, Ethnicity %in% L)

ggplot(wiebe, aes(x = Ethnicity, y= Percentage, fill = wiebe_risk)) +
  geom_bar(stat = "identity")+
  ggtitle("Wiebe et Al. 
DR|DQ Single Molecule Eplet Mismatch Count Risk Categories")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 22))+
  theme(legend.title = element_text(face= "bold", size=20))+
  theme(legend.text = element_text(face= "bold", size=18))+
  theme(axis.title = element_text(face="bold", size= 20))+
  theme(axis.text = element_text(face="bold", size= 18))
ggsave(filename="Wiebe_etAl_DRDQ_ethnicity_distribution_stacked.pdf", width = 10, height = 10, units="in", path = "./")







#hlaR


pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","DON_RACE","DR","DQ","hlaR_risk")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, hlaR_risk == "low")
  ln <- nrow(Low)
  pln <- ln/N*100
  med <- filter(pops, hlaR_risk == "interm")
  mn <- nrow(med)
  pmn <- mn/N*100
  High <- filter(pops, hlaR_risk == "high")
  hn <- nrow(High)
  phn <- hn/N*100
  
  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }
  
  
  Percentage <- c(pln,pmn,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("Low", "Medium", "High")
  pop_row_df$hlaR_risk <- risks
  pop_table<- rbind(pop_table, pop_row_df)
  
  Counts <- c(ln,mn,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("Low", "Medium", "High", "Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_hlaR_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_hlaR_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_hlaR_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)

#pop_table$hlaR_risk <- factor(pop_table$hlaR_risk, levels = c("High", "Medium", "Low"))

ggplot(data = pop_table, aes( x = hlaR_risk, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("hlaR DR|DQ Single Molecule Eplet Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="hlaR_DRDQ_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = hlaR_risk)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("hlaR DR|DQ Single Molecule Eplet Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="hlaR_DRDQ_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")






#I1
pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","DON_RACE","DR","DQ","I1_risk")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, I1_risk == "low")
  ln <- nrow(Low)
  pln <- ln/N*100
  High <- filter(pops, I1_risk == "high")
  hn <- nrow(High)
  phn <- hn/N*100
  
  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }
  
  
  Percentage <- c(pln,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("Low","High")
  pop_row_df$I1_risk <- risks
  pop_table<- rbind(pop_table, pop_row_df)
  
  Counts <- c(ln,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("Low", "High", "Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_Disparity_I1_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_I1_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_I1_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)


ggplot(data = pop_table, aes( x = I1_risk, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("FIBERS I1 AA-MM  Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="I1_risk_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = I1_risk)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("FIBERS I1 AA-MM  Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="I1_risk_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")




#I9
pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","DON_RACE","DR","DQ","I9_risk")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, I9_risk == "low")
  ln <- nrow(Low)
  pln <- ln/N*100
  High <- filter(pops, I9_risk == "high")
  hn <- nrow(High)
  phn <- hn/N*100
  
  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }
  
  
  Percentage <- c(pln,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("Low","High")
  pop_row_df$I9_risk <- risks
  pop_table<- rbind(pop_table, pop_row_df)
  
  Counts <- c(ln,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("Low", "High", "Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_Disparity_I9_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_I9_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_I9_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)


ggplot(data = pop_table, aes( x = I9_risk, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("FIBERS I9 AA-MM  Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="I9_risk_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = I9_risk)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("FIBERS I9 AA-MM  Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="I9_risk_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")




#g1g2
pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","DON_RACE","DR","DQ","g1g2_risk")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, g1g2_risk == "low")
  ln <- nrow(Low)
  pln <- ln/N*100
  med <- filter(pops, g1g2_risk == "interm")
  mn <- nrow(med)
  pmn <- mn/N*100
  High <- filter(pops, g1g2_risk == "high")
  hn <- nrow(High)
  phn <- hn/N*100
  
  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }
  
  
  Percentage <- c(pln,pmn,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("Low", "Medium", "High")
  pop_row_df$g1g2_risk <- risks
  pop_table<- rbind(pop_table, pop_row_df)
  
  Counts <- c(ln,mn,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("Low", "Medium", "High", "Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_Disparity_g1g2_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_g1g2_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_g1g2_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)


ggplot(data = pop_table, aes( x = g1g2_risk, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("DQ G1G2 Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="g1g2_risk_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = g1g2_risk)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("DQ G1G2 Risk Categories
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="g1g2_risk_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")





#REC_DR_MM_EQUIV_CUR
pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","REC_DR_MM_EQUIV_CUR")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, REC_DR_MM_EQUIV_CUR == 0)
  ln <- nrow(Low)
  pln <- ln/N*100
  med <- filter(pops, REC_DR_MM_EQUIV_CUR == 1)
  mn <- nrow(med)
  pmn <- mn/N*100
  High <- filter(pops, REC_DR_MM_EQUIV_CUR == 2)
  hn <- nrow(High)
  phn <- hn/N*100
  
  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }
  
  
  Percentage <- c(pln,pmn,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("AngMM_0", "AngMM_1", "AngMM_2")
  pop_row_df$REC_DR_MM_EQUIV_CUR <- risks
  pop_table<- rbind(pop_table, pop_row_df)
  
  Counts <- c(ln,mn,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("AngMM_0", "AngMM_1", "AngMM_2","Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_Disparity_DRang_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_DRang_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_DRang_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)


ggplot(data = pop_table, aes( x = REC_DR_MM_EQUIV_CUR, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("Recipient DR Antigen Mismatch 
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="REC_DR_MM_EQUIV_CUR_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = REC_DR_MM_EQUIV_CUR)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("Recipient DR Antigen Mismatch
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="REC_DR_MM_EQUIV_CUR_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")



#REC_DQ_MM_EQUIV_CUR
pop_stats <- molecule%>% dplyr::select("PX_ID","CAN_RACE","REC_DQ_MM_EQUIV_CUR")

pop_table= data_frame()
counts_table = data_frame()
stacked_table= data_frame()

for (pop in ethnicity){
  pops <- filter(pop_stats, CAN_RACE== pop)
  N <- nrow(pops)
  Low <- filter(pops, REC_DQ_MM_EQUIV_CUR == 0)
  ln <- nrow(Low)
  pln <- ln/N*100
  med <- filter(pops, REC_DQ_MM_EQUIV_CUR == 1)
  mn <- nrow(med)
  pmn <- mn/N*100
  High <- filter(pops, REC_DQ_MM_EQUIV_CUR == 2)
  hn <- nrow(High)
  phn <- hn/N*100
  
  if (pop == "16: Black or African American"){
    name <- paste0("AFA")
  }
  if (pop == "64: Asian"){
    name <- paste0("ASN")
  }
  if (pop == "8: White"){
    name <- paste0("CAU")
  }
  if (pop == "2000: Hispanic/Latino"){
    name <- paste0("HIS")
  }
  if (pop == "128: Native Hawaiian or Other Pacific Islander"){
    name <- paste0("HPI")
  }
  if (pop == "Multi-Racial"){
    name <- paste0("MLT")
  }
  if (pop == "32: American Indian or Alaska Native"){
    name <- paste0("NAM")
  }
  
  
  Percentage <- c(pln,pmn,phn)
  pop_row_df <- as.data.frame(Percentage)
  pop_row_df$Ethnicity <- name
  risks <- c("AngMM_0", "AngMM_1", "AngMM_2")
  pop_row_df$REC_DQ_MM_EQUIV_CUR <- risks
  pop_table<- rbind(pop_table, pop_row_df)
  
  Counts <- c(ln,mn,hn,N)
  count_row_df <- as.data.frame(Counts)
  count_row_df$Ethnicity <- name
  risky <- c("AngMM_0", "AngMM_1", "AngMM_2","Total")
  count_row_df$Category <- risky
  counts_table<- rbind(counts_table, count_row_df)
}

file_name <- "Ethnic_Disparity_DQang_stats.csv"
write.table(pop_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_DQang_stacked.csv"
write.table(stacked_table,file= file_name, sep = ",", row.names = F)

file_name <- "Ethnic_Disparity_DQang_counts.csv"
write.table(counts_table,file= file_name, sep = ",", row.names = F)


ggplot(data = pop_table, aes( x = REC_DQ_MM_EQUIV_CUR, y = Percentage, fill = Ethnicity ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge' )+
  theme_bw()+
  xlab("Risk Category")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("Recipient DQ Antigen Mismatch 
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))
ggsave(filename="REC_DQ_MM_EQUIV_CUR_ethnicity_distribution_clustered.pdf", width = 10, height = 10, units="in", path = "./")


ggplot(pop_table, aes(x = Ethnicity, y= Percentage, fill = REC_DQ_MM_EQUIV_CUR)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  xlab("Recipient Ethnicity")+
  ylab("Percent of Transplant Recipients")+
  ggtitle("Recipient DQ Antigen Mismatch
by Recipient Ethnic Population")+
  theme(plot.title = element_text(hjust = 0.5, face="bold",size = 16))+
  theme(legend.title = element_text(face= "bold", size=14))+
  theme(legend.text = element_text(face= "bold", size=12))+
  theme(axis.title = element_text(face="bold", size= 14))+
  theme(axis.text = element_text(face="bold", size= 12))+
  coord_flip()
ggsave(filename="REC_DQ_MM_EQUIV_CUR_ethnicity_distribution_stacked_risk.pdf", width = 10, height = 10, units="in", path = "./")







