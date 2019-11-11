library(tidyverse)
library(janitor)
library(here)
library(jntools) 
library(ape)
library(phangorn)
library(ips)
library(gbfetch)

# Process raw specimen data ----

# This includes collection data of all J. Nitta collections (and more).
# Filter this down to just the relevant columns for sporophytes from Nectandra.
specimens_raw <- read_csv(here("data_raw/specimens.csv")) %>%
  clean_names

nectandra_specimens <-
specimens_raw %>%
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
dna_raw <- read_csv(here("data_raw/DNA_accessions.csv")) %>%
  clean_names

nectandra_dna <- dna_raw %>%
  select(genomic_id, specimen_id) %>%
  filter(specimen_id %in% nectandra_specimens$specimen_id)

write_csv(nectandra_dna, "data/nectandra_DNA_accessions.csv")

# Process DNA alignment ----

# Read in all sporophyte fasta sequences.
# These are the exact same sequences that will be submitted to GenBank
# (primers have been trimmed from ends, and sequences annotated with metadata).
nectandra_rbcL <- read.FASTA("data_raw/geneious_sporos_rbcL.fasta")

# Check for missing taxa from Nectandra list.
nectandra_rbcL_missing_taxa <-
tibble(genomic_id = names(nectandra_rbcL)) %>%
  left_join(select(dna_raw, genomic_id, specimen_id)) %>%
  left_join(
    transmute(
      specimens_raw, 
      specimen_id,
      taxon = paste3(genus, species, infraspecific_name)
      )
  ) %>%
  anti_join(nectandra_specimens, ., by = "taxon") %>%
  pull(taxon) %>%
  unique %>% sort

# Construct query and fetch GenBank sequences for missing taxa
query <- glue::glue("({paste(nectandra_rbcL_missing_taxa, collapse = '[ORIGIN] OR ')}[ORIGIN]) AND rbcl")

gb_seqs <- gbfetch::fetch_sequences(query)

gb_seqs_metadata <- gbfetch::fetch_metadata(query)

# Select the appropriate sequences to use and rename them.
# KM008147 = Pteris altissima
# AY175794 = Trichomanes polypodioides
gb_seqs <- gb_seqs[names(gb_seqs) %in% c("KM008147", "AY175794")]

names(gb_seqs) <-
  tibble(accession = names(gb_seqs)) %>%
  left_join(gb_seqs_metadata) %>%
  mutate(tip_label = paste(species, accession) %>% str_replace_all(" ", "_")) %>%
  pull(tip_label)

# Reassign names to be taxon plus DNA accession number
names(nectandra_rbcL) <-
tibble(genomic_id = names(nectandra_rbcL)) %>%
  left_join(select(dna_raw, genomic_id, specimen_id)) %>%
  left_join(
    transmute(
      specimens_raw, 
      specimen_id,
      taxon = paste3(genus, species, infraspecific_name)
    )) %>%
  assert(is_uniq, genomic_id, specimen_id) %>%
  assert(not_na, everything()) %>%
  mutate(tip_label = paste(taxon, genomic_id) %>% str_replace_all(" ", "_")) %>%
  pull(tip_label)

# Combine new sequences with GenBank sequences
nectandra_rbcL <- c(nectandra_rbcL, gb_seqs)

# Align sequences with MAFFT, trim ends, and remove any empty cells
nectandra_rbcL <- mafft(nectandra_rbcL, exec = "/usr/bin/mafft") %>%
  trimEnds(nrow(.) * 0.5) %>%
  deleteEmptyCells()

# Write out alignment in PHYLIP format for RAxML
write.phyDat(nectandra_rbcL, file = "data/nectandra_rbcL.phy")
