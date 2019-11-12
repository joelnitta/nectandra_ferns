# This script takes various raw data files that include many more specimens
# and data columns than those needed for just this project and filters
# and cleans them. Processed data files are written to "data/".

library(janitor)
library(here)
library(jntools) 
library(ape)
library(phangorn)
library(ips)
library(gbfetch)
library(assertr)
library(tidyverse)

# Process raw specimen data ----
#
# Raw specimen data include collection data of all J. Nitta collections (and more).
# Filter this down to just the relevant columns for sporophytes from Nectandra.
# AND two extra specimens that are not from Nectandra but were used for DNA sequencing.

# Read in raw specimen data
specimens_raw <- read_csv(here("data_raw/specimens.csv")) %>%
  clean_names

# Read in and clean up taxonomic data
sci_names <-
read_csv(here("data_raw/taxonomy.csv")) %>%
  janitor::clean_names() %>%
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
  mutate(scientific_name = stringr::str_trim(scientific_name))

# Format final specimen data
nectandra_specimens <-
specimens_raw %>%
  mutate(coll_num = paste3(collection_number, subcollection_number, sep = "")) %>%
  mutate(specimen = paste3(collector_lastname, coll_num)) %>%
  filter(is_gametophyte == 0) %>%
  mutate(notes = replace_na(notes, "")) %>%
  # Filter out specimens that are missing vouchers
  filter(str_detect(notes, "MISSING", negate = TRUE)) %>%
  # Filter to only Nectandra plus two exceptions
  filter(locality == "Nectandra Cloud Forest Preserve" | specimen %in% c("Nitta 2012", "Nitta 169")) %>% 
  select(specimen_id, specimen, 
         genus, species, contains("infraspecific"), certainty,
         country, locality, site, observations,
         elevation_m, latitude, longitude,
         date_collected) %>%
  rename(specific_epithet = species, elevation = elevation_m) %>%
  mutate(species = paste3(genus, specific_epithet)) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  # Remove one Didymoglossum specimen that was separated out from a single collection and
  # may or may not be a distinct species pending additional study.
  filter(taxon != "Didymoglossum sp") %>%
  left_join(sci_names, by = "taxon") %>%
  assert(not_na, specimen_id, specimen, date_collected, genus, species, country, locality)

write_csv(nectandra_specimens, "data/nectandra_specimens.csv")

# Process raw DNA accession data ----

# This includes DNA accession data of all J. Nitta collections (and more).
# Filter this down to just the relevant columns for sporophytes from Nectandra.
# AND two extra specimens that are not from Nectandra but were used for DNA sequencing.
dna_raw <- read_csv(here("data_raw/DNA_accessions.csv")) %>%
  clean_names

nectandra_dna <- dna_raw %>%
  select(genomic_id, specimen_id) %>%
  filter(specimen_id %in% nectandra_specimens$specimen_id)

write_csv(nectandra_dna, "data/nectandra_DNA_accessions.csv")
