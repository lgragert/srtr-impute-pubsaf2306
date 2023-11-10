if (!require("tidyverse"))
  install.packages("tidyverse")
if (!require("haven"))
  install.packages("haven")
if (!require("labelled"))
  install.packages("labelled")
if (!require("dplyr"))
  install.packages("dplyr")

# rec_histo_decoded <- haven::read_sas("rec_histo.sas7bdat", catalog_file = "formats.sas7bcat")

rec_histo_data_file <- "rec_histo.sas7bdat"
catalog_file <- "formats.sas7bcat"

rec_histo <- haven::read_sas(rec_histo_data_file, catalog_file = "formats.sas7bcat")

rec_histo_decoded <- rec_histo %>% mutate (REC_A1 = to_character(REC_A1))

# TODO - determine if standard catalog file includes DR

# summarize character data
table(rec_histo_decoded$REC_A1)


# output to 