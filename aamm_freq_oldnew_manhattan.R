
if (!require("tidyverse"))
  install.packages("tidyverse")
  install.packages("ggtext")
if (!require('ggplot2'))
  install.packages('ggplot2')
if (!require('ggrepel'))
  install.packages('ggrepel')

# Start by getting the old and new AAMM frequency
setwd("./")

new_data <- read.csv('new_AAMM_frequency.csv')
old_data <- read.csv('old_AAMM_frequency.csv')

# GGplot function found on https://r-graph-gallery.com/101_Manhattan_plot.html
manhattan_plot <- function(sub_typ, which_mm, tolerance=0.4) {
  # Make Loci names L1, L2, L3,... so that you can have it in genomic order
  sub_typ$lgroup <- 'NA'
  sub_typ[is.element(sub_typ$Locus,'A'),"lgroup"] <- 'L1'
  sub_typ[is.element(sub_typ$Locus,'C'),"lgroup"] <- 'L2'
  sub_typ[is.element(sub_typ$Locus,'B'),"lgroup"] <- 'L3'
  sub_typ[is.element(sub_typ$Locus,'DRB345'),"lgroup"] <- 'L4'
  sub_typ[is.element(sub_typ$Locus,'DRB1'),"lgroup"] <- 'L5'
  sub_typ[is.element(sub_typ$Locus,'DQA1'),"lgroup"] <- 'L6'
  sub_typ[is.element(sub_typ$Locus,'DQB1'),"lgroup"] <- 'L7'
  sub_typ[is.element(sub_typ$Locus,'DPA1'),"lgroup"] <- 'L8'
  sub_typ[is.element(sub_typ$Locus,'DPB1'),"lgroup"] <- 'L9'
  
  # This makes it to where the y-axis looks nice and the function is customizable
  if (which_mm == 'Overall') {
    sub_typ$Frequency <- sub_typ$Overall_Frequency
  } else if (which_mm == 'One') {
    sub_typ$Frequency <- sub_typ$X1MM_Frequency
  } else {
    sub_typ$Frequency <- sub_typ$X2MM_Frequency
  }
  # Prepare the data
  don <- sub_typ %>%
    
    # Compute locus size by AA position
    group_by(lgroup) %>%
    summarise(chr_len=max(Position)) %>%
    
    # Calculate cumulative position of each loci
    mutate(tot=cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(sub_typ, ., by=c('lgroup'='lgroup')) %>%
    
    # Add a cumulative position of each AA position
    arrange(lgroup, Position) %>%
    mutate(Loci=Position+tot) %>%
    
    # Add highlight to big variance in TRS Average
    mutate(is_annotate=ifelse(Frequency>tolerance, "yes", "no"))
  
  # Make the x-axes
  axisdf = don %>%
    group_by(Locus) %>%
    summarize(center=(max(Loci) + min(Loci)) / 2)
  
  max_y_val <- 0.7
  
  # Plot the data
  ggplot(don, aes(x=Loci, y=Frequency)) +
    
    # Show all points
    geom_point(aes(color=as.factor(lgroup)), alpha=0.8, size=1.3) +
    scale_color_manual(values=rep(c("blue", "red"), 22 )) +
    
    # custom X axis:
    scale_x_continuous(label=axisdf$Locus, breaks=axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_val + 0.01) ) +     # remove space between plot area and x axis
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=Position), size=2, label.size=NA, max.overlaps=nrow(subset(don, is_annotate=="yes"))) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    ggtitle(paste("Frequency of ", which_mm, " AA Mismatch per Loci")) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Frequency of any MM does not separate between 1 or 2 MM
manhattan_plot(new_data, "Overall", 0.3)
ggsave("new_aamm_frequency_manhattan.jpg", width=18.00, height=4.00)

manhattan_plot(old_data, "Overall", 0.05)
ggsave("old_aamm_frequency_manhattan.jpg",width=18.00, height=4.00)

# Frequency of only 1 MM
manhattan_plot(new_data, "One", 0.25)
ggsave("new_1aamm_frequency_manhattan.jpg", width=18.00, height=4.00)

manhattan_plot(old_data, "One", 0.05)
ggsave("old_1aamm_frequency_manhattan.jpg",width=18.00, height=4.00)

# Frequency of only 2 MM
manhattan_plot(new_data, "Two", 0.08)
ggsave("new_2aamm_frequency_manhattan.jpg", width=18.00, height=4.00)

manhattan_plot(old_data, "Two", 0.05)
ggsave("old_2aamm_frequency_manhattan.jpg",width=18.00, height=4.00)
