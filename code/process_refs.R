# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

# Clean bibliography and write out.
clean_bib(
  "ms/references_raw.bib", 
  strip_fields = c("abstract", "file", "keywords", 
                   "publisher", 
                   "url", "issn",
                   "isbn", "month", "number")) %>%
  write_lines("ms/references.bib")

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

