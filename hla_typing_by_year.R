# hla_typing_by_year.R

library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)

setwd("./")

hla<-read.csv("tx_ki_hla_9loc.csv")

# split transplant date to get year
hla<-tidyr::separate(hla,"REC_TX_DT",c('REC_TX_YEAR','REC_TX_MONTH','REC_TX_DATE'))


# Function that creates 2 plots of typing by year
one_loci_plots <- function(locus_rec, locus_don, locus){
  # Change any empty strings into Missing
  # Collapse the factor levels into a single level
  hla$don_typed <- fct_collapse(locus_don, 
                                Missing=c(""), 
                                group_other=TRUE)
  hla$rec_typed <- fct_collapse(locus_rec,
                                Missing=c(""),
                                group_other=TRUE)
  
  # Anything categorized as not missing will appear as typed in the df
  hla$don_typed <- recode(hla$don_typed, 'Other'='DON_TYPED')
  hla$rec_typed <- recode(hla$rec_typed, 'Other'='REC_TYPED')
  
  # Get the recipient proportion typed each year
  rec_prop <- hla %>%
    group_by(REC_TX_YEAR,rec_typed) %>%
    summarize(n=length(PERS_ID)) %>%
    mutate(PROP_TYPED = n / sum(n))
  
  # Create a recipient plot
  ggplot(rec_prop,aes(x=REC_TX_YEAR, y=PROP_TYPED, fill=rec_typed)) +
    geom_bar(stat='identity') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") +
    labs(y=paste("PROP_",locus,'_TYPED', sep='')) +
    ggtitle (paste("Recipient HLA-",locus," Typing Proportions by Year in SRTR SAF Kidney", sep='')) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste("SRTR_Recip_HLA-",locus,"_Typed.jpg", sep=''), width=6.25, height=4.20)
  
  # Get donor proportion typed each year
  don_prop <- hla %>%
    group_by(REC_TX_YEAR,don_typed) %>%
    summarize(n=length(PERS_ID)) %>%
    mutate(PROP_TYPED = n / sum(n))
  
  # Create the donor plot
  ggplot(don_prop,aes(x=REC_TX_YEAR, y=PROP_TYPED, fill=don_typed)) +
    geom_bar(stat='identity') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") +
    labs(y=paste('PROP_',locus,'_TYPED', sep='')) +
    ggtitle (paste("Donor HLA-",locus," Typing Proportions by Year in SRTR SAF Kidney", sep='')) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste("SRTR_Donor_HLA-",locus,"_Typed.jpg", sep=''), width=6.25, height=4.20)
}

# Create plots for all 9 loci
one_loci_plots(locus_rec=hla$REC_A1, locus_don=hla$DON_A1, locus='A')
one_loci_plots(locus_rec=hla$REC_B1, locus_don=hla$DON_B1, locus='B')
one_loci_plots(locus_rec=hla$REC_CW1, locus_don=hla$DON_C1, locus='C')
one_loci_plots(locus_rec=hla$REC_DR1, locus_don=hla$DON_DR1, locus='DR')
one_loci_plots(locus_rec=hla$REC_DQW1, locus_don=hla$DON_DQ1, locus='DQ')
one_loci_plots(locus_rec=hla$REC_DQA1, locus_don=hla$DON_DQA1, locus='DQA1')
one_loci_plots(locus_rec=hla$REC_DPW1, locus_don=hla$DON_DP1, locus='DP')
one_loci_plots(locus_rec=hla$REC_DPA1, locus_don=hla$DON_DPA1, locus='DPA1')

# TODO - Create functions for each broad/split, but for now we just have DQ.

# broad / split DQ
hla$don_dq_split <- fct_collapse(hla$DON_DQ1,
                                 MISSING=c(''),
                                 BROAD_DQ1=c('1'),
                                 SPLIT=c('2','3','4','5','6','7','8','9'),
                                 TWO_FIELD=c('02:01','02:02','03:01','03:02','03:03','03:19',
                                             '04:01','04:02','05:01','05:02','06:01','06:02','06:03','06:04','06:09')
                                 # group_other=TRUE
)

fct_count(hla$don_dq_split)

hla$rec_dq_split <- fct_collapse(hla$REC_DQW1,
                                 MISSING=c(''),
                                 BROAD_DQ1=c('1'),
                                 SPLIT=c('2','3','4','5','6','7','8','9'),
                                 TWO_FIELD=c('02:01','02:02','03:01','03:02','03:03','03:19',
                                             '04:01','04:02','05:01','05:02','06:01','06:02','06:03','06:04','06:09')
)

fct_count(hla$rec_dq_split)


don_dq_split_prop <- hla %>%
  group_by(REC_TX_YEAR,don_dq_split) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DQ_TYPED = n / sum(n))

rec_dq_split_prop <- hla %>%
  group_by(REC_TX_YEAR,rec_dq_split) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DQ_TYPED = n / sum(n))


ggplot(don_dq_split_prop,aes(x=REC_TX_YEAR, y=PROP_DQ_TYPED, fill=don_dq_split)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  # theme(legend.position = "none") +
  ggtitle ("Donor HLA-DQ Broad Split Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5,size=10))

ggsave("SRTR_Donor_HLA-DQ_Broad_Split_Typed.jpg",width=6.25, height=4.20)


ggplot(rec_dq_split_prop,aes(x=REC_TX_YEAR, y=PROP_DQ_TYPED, fill=rec_dq_split)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  # theme(legend.position = "none") +
  ggtitle ("Recipient HLA-DQ Broad Split Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5,size=10))

ggsave("SRTR_Recip_HLA-DQ_Broad_Split_Typed.jpg",width=6.25, height=4.20)
