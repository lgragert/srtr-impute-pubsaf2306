# hla_typing_by_year.R

library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)

setwd("./")

hla<-read.csv("tx_ki_hla_9loc.csv")

# split transplant date to get year
hla<-tidyr::separate(hla,"REC_TX_DT",c('REC_TX_YEAR','REC_TX_MONTH','REC_TX_DATE'))

# hla<-dplyr::count(hla,"DON_C1",name="C_Count")

# squash DON_C1 factor levels
hla$don_c_typed <- fct_collapse(hla$DON_C1,
                         Missing=c(""),
                         group_other=TRUE
                         )

hla$rec_c_typed <- fct_collapse(hla$REC_CW1,
                                 Missing=c(""),
                                 group_other=TRUE
                                 )

hla$don_dq_typed <- fct_collapse(hla$DON_DQ1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$rec_dq_typed <- fct_collapse(hla$REC_DQW1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$don_dp_typed <- fct_collapse(hla$DON_DP1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$rec_dp_typed <- fct_collapse(hla$REC_DPW1,
                                 Missing=c(""),
                                 group_other=TRUE
)

# Add dpa1 and dqa1
hla$don_dqa1_typed <- fct_collapse(hla$DON_DQA1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$rec_dqa1_typed <- fct_collapse(hla$REC_DQA1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$don_dpa1_typed <- fct_collapse(hla$DON_DPA1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$rec_dpa1_typed <- fct_collapse(hla$REC_DPA1,
                                 Missing=c(""),
                                 group_other=TRUE
)

hla$don_c_typed <- recode(hla$don_c_typed,'Other'="DON_C_TYPED")
hla$rec_c_typed <- recode(hla$rec_c_typed,'Other'="REC_C_TYPED")
hla$don_dq_typed <- recode(hla$don_dq_typed,'Other'="DON_DQ_TYPED")
hla$rec_dq_typed <- recode(hla$rec_dq_typed,'Other'="REC_DQ_TYPED")
hla$don_dp_typed <- recode(hla$don_dp_typed,'Other'="DON_DP_TYPED")
hla$rec_dp_typed <- recode(hla$rec_dp_typed,'Other'="REC_DP_TYPED")

hla$don_dqa1_typed <- recode(hla$don_dqa1_typed,'Other'="DON_DQA1_TYPED")
hla$rec_dqa1_typed <- recode(hla$rec_dqa1_typed,'Other'="REC_DQA1_TYPED")
hla$don_dpa1_typed <- recode(hla$don_dpa1_typed,'Other'="DON_DPA1_TYPED")
hla$rec_dpa1_typed <- recode(hla$rec_dpa1_typed,'Other'="REC_DPA1_TYPED")


fct_count(hla$DON_C1)

# recip C

rec_c_prop <- hla %>%
  group_by(REC_TX_YEAR,rec_c_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_C_TYPED = n / sum(n))

ggplot(rec_c_prop,aes(x=REC_TX_YEAR, y=PROP_C_TYPED, fill=rec_c_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Recipient HLA-C Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Recip_HLA-C_Typed.jpg", width=6.25, height=4.20)

# recip DQ

rec_dq_prop <- hla %>%
  group_by(REC_TX_YEAR,rec_dq_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DQ_TYPED = n / sum(n))

ggplot(rec_dq_prop,aes(x=REC_TX_YEAR, y=PROP_DQ_TYPED, fill=rec_dq_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Recipient HLA-DQ Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Recip_HLA-DQ_Typed.jpg",width=6.25, height=4.20)

# recip DP

rec_dp_prop <- hla %>%
  group_by(REC_TX_YEAR,rec_dp_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DP_TYPED = n / sum(n))

ggplot(rec_dp_prop,aes(x=REC_TX_YEAR, y=PROP_DP_TYPED, fill=rec_dp_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Recipient HLA-DP Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Recip_HLA-DP_Typed.jpg",width=6.25, height=4.20)

# recip DQA1

rec_dqa1_prop <- hla %>%
  group_by(REC_TX_YEAR,rec_dqa1_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DQA1_TYPED = n / sum(n))

ggplot(rec_dqa1_prop,aes(x=REC_TX_YEAR, y=PROP_DQA1_TYPED, fill=rec_dqa1_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Recipient HLA-DQA1 Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Recip_HLA-DQA1_Typed.jpg",width=6.25, height=4.20)

# recip DPA1

rec_dpa1_prop <- hla %>%
  group_by(REC_TX_YEAR,rec_dpa1_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DPA1_TYPED = n / sum(n))

ggplot(rec_dpa1_prop,aes(x=REC_TX_YEAR, y=PROP_DPA1_TYPED, fill=rec_dpa1_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Recipient HLA-DPA1 Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Recip_HLA-DPA1_Typed.jpg",width=6.25, height=4.20)


# donor C

don_c_prop <- hla %>%
  group_by(REC_TX_YEAR,don_c_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_C_TYPED = n / sum(n))

ggplot(don_c_prop,aes(x=REC_TX_YEAR, y=PROP_C_TYPED, fill=don_c_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Donor HLA-C Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Donor_HLA-C_Typed.jpg",width=6.25, height=4.20)


# donor DQ

don_dq_prop <- hla %>%
  group_by(REC_TX_YEAR,don_dq_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DQ_TYPED = n / sum(n))

ggplot(don_dq_prop,aes(x=REC_TX_YEAR, y=PROP_DQ_TYPED, fill=don_dq_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Donor HLA-DQ Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Donor_HLA-DQ_Typed.jpg",width=6.25, height=4.20)

# donor DP

don_dp_prop <- hla %>%
  group_by(REC_TX_YEAR,don_dp_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DP_TYPED = n / sum(n))

ggplot(don_dp_prop,aes(x=REC_TX_YEAR, y=PROP_DP_TYPED, fill=don_dp_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Donor HLA-DP Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Donor_HLA-DP_Typed.jpg",width=6.25, height=4.20)


# donor DQA1

don_dqa1_prop <- hla %>%
  group_by(REC_TX_YEAR,don_dqa1_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DQA1_TYPED = n / sum(n))

ggplot(don_dqa1_prop,aes(x=REC_TX_YEAR, y=PROP_DQA1_TYPED, fill=don_dqa1_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Donor HLA-DQA1 Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Donor_HLA-DQA1_Typed.jpg",width=6.25, height=4.20)

# donor DPA1

don_dpa1_prop <- hla %>%
  group_by(REC_TX_YEAR,don_dpa1_typed) %>%
  summarize(n=length(PERS_ID)) %>%
  mutate(PROP_DPA1_TYPED = n / sum(n))

ggplot(don_dpa1_prop,aes(x=REC_TX_YEAR, y=PROP_DPA1_TYPED, fill=don_dpa1_typed)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  ggtitle ("Donor HLA-DPA1 Typing Proportions by Year in SRTR SAF Kidney") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("SRTR_Donor_HLA-DPA1_Typed.jpg",width=6.25, height=4.20)


# broad / split DQ

# hla$DON_DQ1 <- recode_factor(hla$DON_DQ1, .missing=NULL)
# hla$REC_DQW1 <- recode_factor(hla$REC_DQW1, .missing=NULL)

# fct_collapse has issues with group_other=TRUE

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

# hla$don_dq_split <- recode(hla$don_dq_split,'Other'="SPLIT")
# hla$rec_dq_split <- recode(hla$rec_dq_split,'Other'="SPLIT")

# hla$rec_dq_broad <- recode(hla$rec_dq_broad,'Other'="SPLIT")

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
