library(tidyverse)
library(janitor)
library(here)
library(jntools) 
library(ape)
library(phangorn)
library(ips)

# Process raw specimen data ----

# This includes collection data of all J. Nitta collections (and more).
# Filter this down to just the relevant columns for sporophytes from Nectandra.
specimens_raw <- read_csv(here("data_raw/specimens.csv")) 

nectandra_specimens <-
specimens_raw %>%
  clean_names() %>%
  filter(is_gametophyte == 0) %>%
  filter(locality == "Nectandra Cloud Forest Preserve") %>% 
  mutate(coll_num = paste3(collection_number, subcollection_number, sep = "")) %>%
  mutate(specimen = paste3(collector_lastname, coll_num)) %>%
  select(specimen_id, specimen, 
         genus, species, contains("infraspecific"), certainty,
         site, observations,
         elevation_m, latitude, longitude,
         date_collected) %>%
  rename(specific_epithet = species, elevation = elevation_m) %>%
  mutate(species = paste3(genus, specific_epithet)) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  # Remove one Didymoglossum specimen that was separated out from a single collection and
  # may or may not be a distinct species pending additional study.
  filter(taxon != "Didymoglossum sp")

write_csv(nectandra_specimens, "data/nectandra_specimens.csv")

# Process raw DNA accession data ----

# This includes DNA accession data of all J. Nitta collections (and more).
# Filter this down to just the relevant columns for sporophytes from Nectandra.
dna_raw <- read_csv(here("data_raw/DNA_accessions.csv")) 

nectandra_dna <- dna_raw %>%
  clean_names %>%
  select(genomic_id, specimen_id) %>%
  filter(specimen_id %in% nectandra_specimens$specimen_id)

write_csv(nectandra_dna, "data/nectandra_DNA_accessions.csv")

# Process DNA alignment ----

# Read in all sporophyte fasta sequences.
rbcL_raw <- read.FASTA("data_raw/geneious_sporos_rbcL.fasta")

# These include some non-Nectandra specimens, so filter those out.
nectandra_rbcL <- rbcL_raw[names(rbcL_raw) %in% nectandra_dna$genomic_id]

# Reassign names to be taxon plus DNA accession number
names(nectandra_rbcL) <-
tibble(genomic_id = names(nectandra_rbcL)) %>%
  left_join(nectandra_dna) %>%
  left_join(select(nectandra_specimens, specimen_id, taxon)) %>%
  mutate(accession = paste(taxon, genomic_id) %>% str_replace_all(" ", "_")) %>%
  pull(accession)

# Align sequences with MAFFT, trim ends, and remove any empty cells
nectandra_rbcL <- mafft(nectandra_rbcL, exec = "/usr/bin/mafft") %>%
  trimEnds(nrow(.) * 0.5) %>%
  deleteEmptyCells()

# Write out alignment in PHYLIP format for RAxML
write.phyDat(nectandra_rbcL, file = "data/nectandra_rbcL.phy")
