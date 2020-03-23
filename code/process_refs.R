# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

make_ref_list(
  rmd_file = "ms/nectandra_pteridos.Rmd", 
  raw_bib = "ms/references_raw.bib",
  final_bib = "ms/references.bib",
  strip_fields = "abstract",
  exclude = c("Thiers2020", "oet2019", "Leon1992"))

# Make some manual fixes to authors in SI bibliography
# (these are institutions, so need double brackets to
# avoid latex thinking they have first and last names)
read_lines("ms/references.bib") %>%
  str_replace(
    "Pteridophyte Phylogeny Group I",
    "\\{Pteridophyte Phylogeny Group I\\}") %>%
  str_replace(
    "R Core Team",
    "\\{R Core Team\\}") %>%
  str_replace(
    "CBOL Plant Working Group",
    "\\{CBOL Plant Working Group\\}") %>%
  write_lines("ms/references.bib")

# To edit the PLoS ONE style, go to
# https://csl.mendeley.com/ and use the visual editor.
# Save the CSL file as 'plos-jhn-2'
# Download to the 'ms' folder:
read_lines("https://csl.mendeley.com/styles/25428611/plos-jhn-2") %>%
  write_lines("ms/plos-one.csl")
