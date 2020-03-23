# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

make_ref_list(
  rmd_file = "ms/nectandra_pteridos.Rmd", 
  raw_bib = "ms/references_raw.bib",
  final_bib = "ms/references.bib",
  strip_fields = "abstract",
  exclude = c("Thiers2020", "oet2019"))

# To edit the PLoS ONE style, go to
# https://csl.mendeley.com/ and use the visual editor.
# Save the CSL file as 'plos-jhn-2'
# Download to the 'ms' folder:
read_lines("https://csl.mendeley.com/styles/25428611/plos-jhn-2") %>%
  write_lines(csl, "ms/plos-one.csl")
