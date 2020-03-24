# This script takes various raw data files that include many more specimens
# and data columns than those needed for just this project and filters
# and cleans them. Processed data files are written to "data/".

library(lubridate)
library(janitor)
library(here)
library(jntools) 
library(ape)
library(phangorn)
library(ips)
library(gbfetch)
library(assertr)
library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")

# Process raw specimen data ----

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
  # Manipulate columns
  mutate(coll_num = paste3(collection_number, subcollection_number, sep = "")) %>%
  mutate(specimen = paste3(collector_lastname, coll_num)) %>%
  mutate(collector = paste(collector_firstname, collector_lastname)) %>%
  rename(specific_epithet = species, elevation = elevation_m, other_collectors = collectors_other) %>%
  mutate(species = paste3(genus, specific_epithet)) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  mutate(notes = replace_na(notes, "")) %>%
  # Format herbaria labels to follow Index Herbariorum (http://sweetgum.nybg.org/science/ih/)
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
  # Format collection date YYYY-MM-DD
  mutate(
    date_collected = date_collected %>%
      str_remove_all("^[:alpha:]+ ") %>% 
      str_remove_all("00:00:00 [:alpha:]+ ") %>%
      mdy() %>%
      as_date(),
    year = year(date_collected),
    month = month(date_collected) %>% str_pad(side = "left", pad = "0", width = 2),
    day = day(date_collected) %>% str_pad(side = "left", pad = "0", width = 2),
    date_collected = paste(year, month, day, sep = "-")
  ) %>%
  # Filter out gametophytes
  filter(is_gametophyte == 0) %>%
  # Filter out specimens that are missing vouchers
  filter(str_detect(notes, "MISSING", negate = TRUE)) %>%
  # Filter to only Nectandra plus two exceptions
  filter(locality == "Nectandra Cloud Forest Preserve" | specimen %in% c("Nitta 2012", "Nitta 169")) %>% 
  # Remove one Didymoglossum specimen that was separated out from a single collection and
  # may or may not be a distinct species pending additional study.
  filter(taxon != "Didymoglossum sp") %>%
  # Add scientific names
  left_join(sci_names, by = "taxon") %>% 
  # Select variables
  select(specimen_id, specimen, 
         genus, specific_epithet, infraspecific_rank, infraspecific_name, certainty,
         species, taxon, scientific_name,
         country, locality, site, observations,
         elevation, latitude, longitude,
         collector, other_collectors,
         herbaria,
         date_collected) %>%
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
  filter(specimen_id %in% nectandra_specimens$specimen_id) %>%
  assert(not_na, genomic_id, specimen_id) %>%
  assert(is_uniq, genomic_id)

write_csv(nectandra_dna, "data/nectandra_DNA_accessions.csv")

# Process Costa Rica pterido richness data ----

# Count number of taxa at Nectandra

  
# From the raw flora list of La Selva (Feb. 2017),
# subset to only pteridophytes and count the number of 
# unique taxa.

la_selva_pteridos <- readxl::read_excel("data_raw/Lista_especies_LS_feb2017.xlsx", skip = 12) %>%
  janitor::clean_names() %>%
  rename(family = familia, genus = genero, specific_epithet = especie, 
         author = autor, notes = historia_taxonomica,
         habit = habito, habit_atr = habito_atributo) %>%
  left_join(taxize::apgFamilies(), by = "family") %>%
  filter(order %in% c(
    "Cyatheales",
    "eupolypod II",
    "Gleicheniales" ,
    "Hymenophyllales",
    "Lycopodiales",
    "Marattiales",
    "Ophioglossales",
    "Polypodiales",
    "Polypodiales-eupolypod I",
    "Polypodiales-eupolypod II"
  )) %>%
  mutate(taxon = jntools::paste3(genus, specific_epithet))

costa_rica_richness_raw <- read_csv("data_raw/costa_rica_richness_raw.csv")

costa_rica_richness <-
  costa_rica_richness_raw %>%
  mutate(richness = as.integer(richness)) %>%
  mutate(richness = case_when(
    name == "La Selva" ~ n_distinct(la_selva_pteridos$taxon),
    name == "Nectandra" ~ n_distinct(nectandra_specimens$species),
    TRUE ~ richness))

write_csv(costa_rica_richness, "data/costa_rica_richness.csv")

# Process PPGI data ----

# Read in original Pteridophyte Phylogeny Group I taxonomy scheme
# (genus-level and higher), add updates, and keep only needed columns.

ppgi_raw <- read_csv("data_raw/ppgi_taxonomy_raw.csv")

ppgi <- ppgi_raw %>%
  # Update with Hiya
  bind_rows(
    tibble(
      class = "Polypodiopsida",
      order = "Polypodiales",
      suborder = "Dennstaedtiineae",
      order_id = "N",
      fam_id = 31,
      family = "Dennstaedtiaceae",
      genus = "Hiya"
    )
  ) %>%
  # Keep only needed columns
  select(class, order, suborder, family, subfamily, genus)

write_csv(ppgi, "data/ppgi_taxonomy.csv")
