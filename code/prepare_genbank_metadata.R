# Script to prepare metadata for submission with GenBank sequences
# via Geneious

library(lubridate)
library(assertr)
library(janitor)
library(jntools)
library(tidyverse)
library(here)
library(glue)
library(stringr)
library(writexl)
library(conflicted)

conflict_prefer("here", "here")
conflict_prefer("filter", "dplyr")

# Read in list of scientific names with authors
sci_names <- read_csv(here("data_raw/taxonomy.csv")) %>% 
  clean_names() %>%
  rename(specific_epithet = species) %>%
  mutate(taxon = paste3(genus, specific_epithet, infrasp_name)) %>%
  mutate(scientific_name = paste3(
    genus,
    specific_epithet,
    author,
    infrasp_rank,
    infrasp_name,
    var_author
  )) %>%
  select(taxon, scientific_name) %>%
  mutate(scientific_name = str_trim(scientific_name)) %>%
  unique() %>%
  assert(is_uniq, taxon)

# Read in list of all DNA accessions
dna_raw <- read_csv(here("data_raw/DNA_accessions.csv")) %>%
  clean_names() %>%
  select(genomic_id, specimen_id) %>%
  assert(is_uniq, genomic_id) %>%
  assert(not_na, genomic_id, specimen_id)

# Read in specimen collection data
specimens_raw <- read_csv(here("data_raw/specimens.csv")) %>%
  clean_names() %>%
  rename(specific_epithet = species) %>%
  mutate(taxon_with_rank = paste3(genus, specific_epithet, infraspecific_rank, infraspecific_name)) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name))

# Combine specimen, DNA accession, and taxonomy data into metadata for GenBank submission
genbank_metadata <- 
  specimens_raw %>%
  mutate(locality = case_when(
    locality == "Palmira" ~ "Palmira, Alajuela Prov.",
    TRUE ~ locality
  )) %>%
  transmute(
    specimen_id, 
    organism = taxon_with_rank,
    taxon,
    country = paste(country, locality, sep = ": "),
    collection_number = paste3(collection_number, subcollection_number, sep = ""),
    collector_lastname,
    herbaria,
    date_collected) %>%
  mutate(
    herbaria = str_remove_all(herbaria, "Nectandra Cloud Forest Preserve") %>%
      str_remove_all("Nectandra") %>%
      str_remove_all("TI") %>%
      str_replace_all("INB", "CR") %>%
      str_replace_all("CR, CR", "CR") %>%
      str_trim %>%
      str_remove_all(", ,") %>%
      str_remove_all(",$")
  ) %>% 
  mutate(voucher = glue("{collector_lastname} {collection_number} ({herbaria})")) %>%
  inner_join(dna_raw) %>%
  left_join(sci_names) %>%
  mutate(
    date_collected = date_collected %>%
      str_remove_all("^[:alpha:]+ ") %>% 
      str_remove_all("00:00:00 [:alpha:]+ ") %>%
      mdy() %>%
      as_date(),
    year = year(date_collected),
    month = month(date_collected, label = TRUE) %>% as.character(),
    day = day(date_collected) %>% str_pad(side = "left", pad = "0", width = 2),
    collection_date = paste(day, month, year, sep = "-")
  ) %>% 
  assert(not_na, genomic_id) %>%
  assert(is_uniq, genomic_id) %>%
  select(g_number = genomic_id, species = organism, country, Voucher = voucher, note = scientific_name, collection_date)

# This will need to be converted to "xls" with excel for use by 
# the Geneious biocode plugin
writexl::write_xlsx(genbank_metadata, here("data_raw/genbank_metadata.xlsx"))
