# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

make_ref_list(
  rmd_file = "ms/nectandra_pteridos.Rmd", 
  raw_bib = "ms/references_raw.bib",
  final_bib = "ms/references.bib",
  strip_fields = "abstract")
