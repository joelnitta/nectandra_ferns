### Set up bibliography and citation style files ###

library(tidyverse)
library(assertr)
library(RefManageR)

# Filter bibtex bibliography ----

# Set path to RMD file
rmd_file <- "ms/nectandra_ferns.Rmd"

# Parse RMD file and extract citation keys
citations <- 
  readr::read_lines(rmd_file) %>% 
  stringr::str_split(" |;") %>% 
  unlist %>% 
  magrittr::extract(., stringr::str_detect(., "@")) %>% 
  stringr::str_remove_all("\\[|\\]|\\)|\\(|\\.$|,|\\{|\\}") %>% 
  magrittr::extract(., stringr::str_detect(., "^@")) %>% 
  stringr::str_remove_all("@") %>% 
  unique %>% sort

# Read in entire reference library exported from Zotero with BetterBibTex
bib <- RefManageR::ReadBib("data_raw/main_library.bib", check = FALSE)

# Convert reference library to dataframe
bib_df <- as.data.frame(bib) %>% rownames_to_column("cite_key")

# Make sure no citations are missing
tibble(cite_key = citations) %>%
  anti_join(bib_df, by = "cite_key") %>%
  verify(nrow(.) == 0, success_fun = success_logical)

# Filter bibliography to only cited references, write it out
bib_filtered <- bib_df %>%
  filter(cite_key %in% citations) %>%
  # RefManageR::as.BibEntry() needs rownames to be citation keys
  column_to_rownames("cite_key") %>%
  RefManageR::as.BibEntry() %>%
  RefManageR::WriteBib(file = "ms/references.bib")

# Download edited PLoS ONE citation style ----

# To edit the PLoS ONE style, go to
# https://csl.mendeley.com/ and use the visual editor.
# Save the CSL file as 'plos-jhn-2'
# Download to the 'ms' folder:
read_lines("https://csl.mendeley.com/styles/25428611/plos-jhn-2") %>%
  write_lines("ms/plos-one.csl")

