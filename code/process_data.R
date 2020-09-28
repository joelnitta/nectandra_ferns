# This script takes various raw data files that include many more specimens
# and data columns than those needed for just this project and filters
# and cleans them. Processed data files are written to "data/".

library(lubridate)
library(janitor)
library(here)
library(jntools) 
library(assertr)
library(tidyverse)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")

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
  mutate(scientific_name = stringr::str_trim(scientific_name)) %>%
  assert(is_uniq, scientific_name, taxon) %>%
  assert(not_na, scientific_name, taxon)

# Format final specimen data
nectandra_specimens <-
  specimens_raw %>%
  # Manipulate columns
  mutate(coll_num = paste3(collection_number, subcollection_number, sep = "")) %>%
  mutate(specimen = paste3(collector_lastname, coll_num)) %>%
  mutate(collector = paste(collector_firstname, collector_lastname)) %>%
  rename(elevation = elevation_m, other_collectors = collectors_other) %>%
  mutate(species = paste3(genus, specific_epithet)) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  mutate(notes = replace_na(notes, "")) %>%
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
  mutate(is_gametophyte = replace_na(is_gametophyte, 0)) %>% # assume NA means sporophyte
  assert(assertr::in_set(c(0,1)), is_gametophyte) %>%
  filter(is_gametophyte == 0) %>%
  # Filter out specimens that are missing vouchers
  filter(str_detect(notes, "MISSING", negate = TRUE)) %>%
  # Filter to only Nectandra
  filter(locality == "Nectandra Cloud Forest Preserve") %>% 
  # Remove one Didymoglossum specimen that was separated out from a single collection and
  # may or may not be a distinct species pending additional study.
  filter(taxon != "Didymoglossum sp") %>%
  # Remove a non-pteridophyte collection
  filter(taxon != "Clethra mexicana") %>%
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
  # Run final checks
  assert(not_na, specimen_id, specimen, date_collected, genus, species, country, locality) %>%
  assert(is_uniq, specimen_id, specimen) %>%
  assert(in_set(ppgi$genus), genus)

write_csv(nectandra_specimens, "data/nectandra_specimens.csv")

# Process raw DNA accession data ----

# This includes DNA accession data of all J. Nitta collections (and more).
# Filter this down to just the relevant columns for sporophytes from Nectandra.
dna_raw <- read_csv(here("data_raw/DNA_accessions.csv")) %>%
  clean_names

nectandra_dna <- dna_raw %>%
  select(genomic_id, specimen_id) %>%
  filter(specimen_id %in% nectandra_specimens$specimen_id) %>%
  assert(not_na, genomic_id, specimen_id) %>%
  assert(is_uniq, genomic_id)

write_csv(nectandra_dna, "data/nectandra_DNA_accessions.csv")

# Process Costa Rica pterido richness data ----

# Count number of pteridophyte taxa at La Selva
# (their list includes hybrids but not varieties)

# Read in raw flora list of La Selva (Feb. 2017)
la_selva_pteridos_raw <- readxl::read_excel("data_raw/Lista_especies_LS_feb2017.xlsx", skip = 12)

# First, filter to only pteridophytes  by using APG orders.
la_selva_pteridos <- 
  la_selva_pteridos_raw %>%
  clean_names() %>%
  rename(family = familia, genus = genero, specific_epithet = especie, 
         author = autor, notes = historia_taxonomica,
         habit = habito, habit_atr = habito_atributo) %>%
  mutate(taxon = jntools::paste3(genus, specific_epithet)) %>%
  # Join to APG order, then filter to only pteridophyte orders
  left_join(taxize::apgFamilies(), by = "family") %>%
  assert(not_na, order) %>%
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
  ))

# As a double-check that I've filtered to only and all pteridophytes at La Selva,
# do an alternate filtering method by filtering to only families in PPGI
la_selva_pteridos_alt <-
  la_selva_pteridos_raw %>%
  clean_names() %>%
  rename(family = familia, genus = genero, specific_epithet = especie, 
         author = autor, notes = historia_taxonomica,
         habit = habito, habit_atr = habito_atributo) %>%
  mutate(taxon = jntools::paste3(genus, specific_epithet)) %>%
  filter(family %in% ppgi$family) %>%
  # Check that the numbers match
  verify(n_distinct(taxon) == n_distinct(la_selva_pteridos$taxon))

# Read in raw richness numbers for other areas from literature
costa_rica_richness_raw <- read_csv("data_raw/costa_rica_richness_raw.csv")

# Exclude non-native species from Nectandra number
nectandra_specimens_native_only <- filter(nectandra_specimens, taxon != "Macrothelypteris torresiana")

# Fill in species richness at Nectandra and La Selva
costa_rica_richness <-
  costa_rica_richness_raw %>%
  mutate(richness = as.integer(richness)) %>%
  mutate(richness = case_when(
    name == "La Selva" ~ n_distinct(la_selva_pteridos_alt$taxon),
    name == "Nectandra" ~ n_distinct(nectandra_specimens_native_only$species),
    TRUE ~ richness)) %>%
  mutate(richness_per_ha = richness / area_ha) %>%
  select(name, full_name, 
         min_el_m, max_el_m, area_ha, 
         richness,  richness_per_ha, 
         holdridge_type, citation, citation_number,
         latitude, longitude)

# Write out final richness data
write_csv(costa_rica_richness, "data/costa_rica_richness.csv")
