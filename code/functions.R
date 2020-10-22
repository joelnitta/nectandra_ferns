# Data wrangling ----

#' Unzip Nitta et al 2017 Ecol Mono data zip file and 
#' extract needed data files.
#' 
#' The dryad data zip file should be downloaded from 
#' https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g
#' (click on "Download dataset")
#'
#' @param dryad_zip_file Path to the data zip file downloaded from Dryad.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return A dataframe describing the unzipped data file.
#' The unzipped data file (rbcL_clean_sporos.fasta: rbcL sequences of 
#' sporophytes from Moorea) will be written to `unzip_path`.
#'
unzip_nitta_2017 <- function (dryad_zip_file, exdir, ...) {
  
  # Dryad data are nested: a zip file inside a zip file.
  # Unzip the first one in a temporary directory, then the one we want into data_raw.
  temp_dir <- tempdir()
  unzip(dryad_zip_file, "data_and_scripts.zip", exdir = temp_dir)
  unzip(fs::path(temp_dir, "data_and_scripts.zip"), files = "data_and_scripts/shared_data/rbcL_clean_sporos.fasta", exdir = exdir, junkpaths = TRUE)
  
  unzipped_files <- unzip(fs::path(temp_dir, "data_and_scripts.zip"), files = "data_and_scripts/shared_data/rbcL_clean_sporos.fasta", list = TRUE)
  
  # Cleanup
  fs::file_delete(fs::path(temp_dir, "data_and_scripts.zip"))
  
  return(unzipped_files)
  
}

#' Unzip Ebihara and Nitta 2017 Ecol Mono data zip file and 
#' extract needed data files.
#' 
#' The dryad data zip file should be downloaded from 
#' https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g
#' (click on "Download dataset")
#'
#' @param dryad_zip_file Path to the data zip file downloaded from Dryad.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return A dataframe of the unzipped files. Externally, these files will be written to `exdir`.
#'
unzip_ebihara_2019 <- function (dryad_zip_file, exdir, ...) {
  
  # Unzip only the needed files
  unzip(dryad_zip_file, files = c("rbcl_mrbayes.nex", "FernGreenListV1.01E.xls", "ESM1.csv"), exdir = exdir, overwrite = TRUE)
 
  # Return dataframe of unzipped files
  unzip(dryad_zip_file, files = c("rbcl_mrbayes.nex", "FernGreenListV1.01E.xls","ESM1.csv"), list = TRUE)
  
}

#' Tidy taxonomic data of pteridophytes of Japan
#'
#' @param Dataframe with data of pteridophytes of Japan (Japan Green list)
#'
tidy_japan_names <- function (data) {
  data %>%
  select(taxon_id = ID20160331, scientific_name = `GreenList Name`,
         endemic = Endemism, conservation_status = RL2012) %>%
    mutate(taxon_id = as.character(taxon_id)) %>%
    select(taxon_id, scientific_name)
}

#' Rename taxa in rbcL alignment of pteridophytes of Japan
#' 
#' Original names are formatted as a series of numbers 
#' separated by underscore, e.g., 601_1. 
#' The first number (601) is a family code, 
#' and the second (1) is the taxon code.
#' 
#' This renames them to human-readable taxon names.
#' 
#' Uses parse_names_batch, which requires GNparser to be installed
#' and on PATH
#'
#' @param japan_rbcL rbcL alignment with names as codes
#' @param japan_taxa dataframe matching codes to species name
#'
rename_japan_rbcL <- function (japan_rbcL, japan_taxa) {
  
  # Make a table mapping dna taxon codes to scientific names
  japan_names_table <-
    tibble(taxon_id = names(japan_rbcL) %>%
             str_split("_") %>%
             map_chr(2)
    ) %>%
    assert(is_uniq, taxon_id) %>%
    assert(not_na, taxon_id) %>%
    left_join(japan_taxa, by = "taxon_id") %>%
    # Convert full name with authors to only taxon, no authors
    assert(is_uniq, scientific_name) %>%
    assert(not_na, scientific_name) %>%
    left_join(
      taxastand::parse_names_batch(.$scientific_name) %>% 
        select(scientific_name = query, taxon),
      by = "scientific_name"
    ) %>%
    assert(is_uniq, taxon) %>%
    assert(not_na, taxon) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
    # Add JA so we know where it came from
    mutate(taxon = paste0(taxon, "_JA"))
  
  # Rename DNA sequences with taxon names
  names(japan_rbcL) <- japan_names_table$taxon
  
  japan_rbcL
}

#' Rename taxa in rbcL alignment of pteridophytes of Japan,
#' and filter to only sexual diploid taxa
#' 
#' Original names are formatted as a series of numbers 
#' separated by underscore, e.g., 601_1. 
#' The first number (601) is a family code, 
#' and the second (1) is the taxon code.
#' 
#' This renames them to human-readable taxon names, and adds
#' "JAsexdip" to each name.
#' 
#' Uses parse_names_batch, which requires GNparser to be installed
#' and on PATH
#'
#' @param japan_rbcL rbcL alignment with names as codes
#' @param japan_taxa dataframe matching codes to species name
#' @param japan_repro_data dataframe with columns including `taxon_id` and `sexual_diploid`
#'
rename_japan_rbcL_sexdip <- function (japan_rbcL, japan_taxa, japan_repro_data) {
  
  # Make a table mapping dna taxon codes to scientific names
  # only 
  japan_sexdip_names_table <-
    tibble(
      original_name = names(japan_rbcL),
      taxon_id = names(japan_rbcL) %>%
        str_split("_") %>%
        map_chr(2)
    ) %>%
    assert(not_na, taxon_id, original_name) %>%
    assert(is_uniq, taxon_id, original_name) %>%
    left_join(japan_taxa, by = "taxon_id") %>%
    # Convert full name with authors to only taxon, no authors
    assert(not_na, scientific_name) %>%
    assert(is_uniq, scientific_name) %>%
    left_join(
      taxastand::parse_names_batch(.$scientific_name) %>% 
        select(scientific_name = query, taxon),
      by = "scientific_name"
    ) %>%
    assert(not_na, taxon) %>%
    assert(is_uniq, taxon) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
    left_join(select(japan_repro_data, taxon_id, sexual_diploid), by = "taxon_id") %>%
    # Filter to only sexual diloids
    mutate(sexual_diploid = replace_na(sexual_diploid, 0)) %>%
    filter(sexual_diploid == 1) %>%
    # Add JA so we know where it came from
    mutate(taxon = paste0(taxon, "_JAsexdip"))
  
  # Subset sequences to only sexual diploids
  japan_rbcL <- japan_rbcL[japan_sexdip_names_table$original_name]
  
  # Rename DNA sequences with taxon names
  names(japan_rbcL) <- japan_sexdip_names_table$taxon
  
  japan_rbcL
}

#' Download altitude data for a country
#'
#' To see list of 3-letter ISO codes, run
#' raster::getData('ISO3')
#'
#' @param country Country name specified by 3-letter ISO code.
#' Default = "CRI" (Costa Rica)
#'
#' @return Dataframe. Columns:
#' - lat: latitidue
#' - lon: longitude
#' - alt: altitude (in m)
#'
download_country_el <- function(country = "CRI") {
  
  dem.raster <- raster::getData("alt", country = country, path = tempdir())
  dem.m  <- raster::rasterToPoints(dem.raster)
  dem.df <- data.frame(dem.m)
  colnames(dem.df) = c("lon", "lat", "alt")
  
  tibble::as_tibble(dem.df)
  
}

#' Load genbank accession numbers received from NCBI after
#' submitting sequences
#'
#' Also adds one more accession for a sequence that
#' was submitted separately (pacific_ferns_rbcL.sqn Nitta2237 MT657442).
#'
#' @param file Path to text file with accession numbers
#'
#' @return Dataframe
#' 
load_genbank_accs <- function (file) {
  
  # Read in seqids.txt received from GenBank after submitting sequences
  read_tsv(file, col_names = c("seq_name", "genbank_accession")) %>%
    # Parse genomicID field
    separate(seq_name, c("file", "genomic_id"), sep = " ") %>%
    select(genomic_id, genbank_accession) %>%
    # Add one more accession that was submitted separately:
    # JNG4254 (Nitta 2237 Amauropelta atrovirens)
    # pacific_ferns_rbcL.sqn Nitta2237       	MT657442
    bind_rows(
      tibble(
        genomic_id = "JNG4254",
        genbank_accession = "MT657442"
      )) %>%
    assert(is_uniq, genomic_id, genbank_accession) %>%
    assert(not_na, genomic_id, genbank_accession)
  
}

# Checklist ----

#' Make a species checklist of pteridophytes at Nectandra
#'
#' @param specimens Specimens data, one line per specimen
#' @param taxonomy Higher level taxonomy data (genus and above)
#'
#' @return Tibble; species checklist of pteridophytes at Nectandra,
#' with one row per taxon
#' 
make_checklist <- function (specimens, taxonomy) {
  
  specimens %>%
    assert(not_na, genus) %>%
    left_join(taxonomy, by = "genus") %>%
    assert(not_na, family) %>%
    mutate(habit_simple = case_when(
      str_detect(observations, regex("epipet", ignore_case = TRUE)) ~ "epipetric",
      str_detect(observations, regex("on rocks", ignore_case = TRUE)) ~ "epipetric",
      str_detect(observations, regex("epiphy", ignore_case = TRUE)) ~ "epiphytic",
      str_detect(observations, regex("growing in bryophyte", ignore_case = TRUE)) ~ "epiphytic",
      str_detect(observations, regex("terr", ignore_case = TRUE)) ~ "terrestrial",
      str_detect(observations, regex("climbing", ignore_case = TRUE)) ~ "climbing",
      str_detect(observations, regex("clambering", ignore_case = TRUE)) ~ "clambering"
    )) %>%
    mutate(scientific_name = case_when(
      specific_epithet == "sp1" ~ paste(genus, specific_epithet),
      TRUE ~ scientific_name
    )) %>%
    group_by(class, family, scientific_name) %>%
    summarize(
      voucher = paste(specimen, collapse = ", "),
      habit = habit_simple %>% 
        magrittr::extract(!is.na(.)) %>% 
        unique %>% sort %>% paste(collapse = ", ") %>%
        stringr::str_to_sentence()
    ) %>%
    ungroup
  
}

# Collection curve ----

#' Estimate species richness using days of sampling as
#' the sampling unit
#'
#' @param specimens Specimens dataframe
#' @param endpoint Endpoint of extrapolation
#'
#' @return Tibble; output of iNEXT()
#' 
estimate_richness_by_date <- function (specimens, endpoint = 150) {

  specimens %>%
    assert(not_na, date_collected) %>%
    # Convert to presence-absence matrix of species x dates
    select(taxon, date_collected) %>%
    unique %>%
    mutate(abun = 1) %>%
    pivot_wider(names_from = date_collected, values_from = abun, 
                values_fill = list(abun = 0)) %>%
    column_to_rownames("taxon") %>% 
    as.matrix %>%
    # Run iNEXT on incidence data
    iNEXT(datatype = "incidence_raw", endpoint = endpoint)

  }

# Barcode analysis ----

#' Make a list of "missing" Nectandra taxa
#' 
#' "Missing" taxa are those that could not be successfully sequenced
#'
#' @param nectandra_rbcL_raw Unaligned rbcL sequences of pteridophytes from Nectandra
#' @param nectandra_dna DNA accession numbers linking accession to specimen ID
#' @param nectandra_specimens Specimen collection data including species names
#' 
#' @return Character vector
make_missing_taxa_list <- function (nectandra_rbcL_raw, nectandra_dna, nectandra_specimens) {
  
  # Check for missing taxa from Nectandra list.
  tibble(genomic_id = names(nectandra_rbcL_raw)) %>%
    left_join(nectandra_dna, by = "genomic_id") %>%
    left_join(select(nectandra_specimens, specimen_id, taxon), by = "specimen_id") %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, genomic_id) %>%
    anti_join(nectandra_specimens, ., by = "taxon") %>%
    pull(taxon) %>%
    unique %>% sort
  
}

#' Tidy GenBank metadata
#'
#' GenBank metadata comes with a column "subtype" containing miscellaneous
#' data separated by '|'. The names of these data are in the "subname" column.
#' This function tidies these two columns, i.e., converts the data in the
#' "subtype" column so each value gets its own column.
#' 
#' @param data Dataframe; GenBank metadata obtained with
#' gbfetch::fetch_metadata.
#'
#' @return Dataframe.
#' 
#' @examples
#' raw_meta <- fetch_metadata("rbcl[Gene] AND Crepidomanes[ORGN]")
#' # Raw metadata still contains untidy data in the "subtype" and 
#' # "subname" columns
#' raw_meta
#' # Tidy!
#' tidy_genbank_metadata(raw_meta)

tidy_genbank_metadata <- function(data) {
  
  # Check assumptions: must have subtype and subname columns present
  assertthat::assert_that("subtype" %in% colnames(data))
  assertthat::assert_that("subname" %in% colnames(data))
  
  # Define helper function that works on one row at a time
  tidy_genbank_meta_single <- function(data) {
    
    # Early return if no data to parse in subtype
    if(data$subtype == "" | is.na(data$subtype)) return (data %>% dplyr::select(-subname, -subtype))
    if(data$subname == "" | is.na(data$subname)) return (data %>% dplyr::select(-subname, -subtype))
    
    # Split "subtype" into multiple columns
    # Use janitor::make_clean_names to de-duplicate names
    sub_cols <- data %>% 
      dplyr::pull(subtype) %>% 
      stringr::str_split("\\|") %>% 
      magrittr::extract2(1) %>% 
      janitor::make_clean_names()
    
    sub_data <- data %>% 
      dplyr::select(subname) %>% 
      tidyr::separate(subname, into = sub_cols, sep = "\\|")
    
    data %>% 
      dplyr::select(-subname, -subtype) %>% 
      dplyr::bind_cols(sub_data)
  }
  
  # Apply the function row-wise and combine back into a dataframe
  transpose(data) %>%
    purrr::map(as_tibble) %>%
    purrr::map_df(tidy_genbank_meta_single)
  
}

#' Download genbank sequences for "missing" Nectandra taxa
#' 
#' "Missing" taxa are those that could not be successfully sequenced
#'
#' @param nectandra_rbcL_missing_taxa Character vector of missing taxa
#' @param acc_keep Character vector of GenBank accessions to use
#' @return List of class "DNAbin"
#' 
fetch_gb_seqs <- function (nectandra_rbcL_missing_taxa, acc_keep) {
  
  # Construct query and fetch GenBank sequences for missing taxa
  query <- glue::glue("({paste(nectandra_rbcL_missing_taxa, collapse = '[ORIGIN] OR ')}[ORIGIN]) AND rbcl")
  
  gb_seqs <- gbfetch::fetch_sequences(query)
  
  gb_seqs_metadata <- gbfetch::fetch_metadata(query)
  
  gb_seqs_metadata <- gb_seqs_metadata %>%
    filter(accession %in% acc_keep) %>%
    tidy_genbank_metadata %>%
    rename(specimen = specimen_voucher) %>%
    # Manually add data where missing
    mutate(specimen = str_remove_all(specimen, "\\(UC\\)|M\\.|\\(Duke\\)|J\\.-Y\\.") %>% str_trim("both")) %>% 
    mutate(country = case_when(
      species == "Radiovittaria remota" ~ "Costa Rica",
      TRUE ~ country
    )) %>%
    mutate(specimen = case_when(
      species == "Radiovittaria remota" ~ "Moran 3180a",
      TRUE ~ specimen
    )) %>%
    mutate(specimen = case_when(
      species != "Pteris altissima" ~ paste(specimen, country),
      TRUE ~ specimen
    )) %>%
    select(taxon = species, specimen, genbank_accession = accession) %>%
    mutate(tip_label = paste(taxon, specimen) %>% str_replace_all(" ", "_"))
  
  # Select the appropriate sequences to use and rename them as species_accession
  gb_seqs <- gb_seqs[names(gb_seqs) %in% acc_keep]
  
  names(gb_seqs) <-
    tibble(genbank_accession = names(gb_seqs)) %>%
    left_join(gb_seqs_metadata, by = "genbank_accession") %>%
    pull(tip_label)
  
  list(
    seqs = gb_seqs,
    metadata = gb_seqs_metadata)
  
}

#' Rename Nectandra rbcL sequences
#' 
#' Renames newly sequenced samples to "Species_Voucher"
#'
#' @param nectandra_rbcL_raw Unaligned rbcL sequences of pteridophytes from Nectandra
#' @param nectandra_dna DNA accession numbers linking accession to specimen ID
#' @param nectandra_specimens Specimen collection data including species names
#'
#' @return List of class DNAbin
#' 
rename_nectandra_rbcL <- function (nectandra_rbcL_raw, nectandra_dna, nectandra_specimens) {
  
  # Reassign names to be taxon plus DNA accession number
  names(nectandra_rbcL_raw) <-
    tibble(genomic_id = names(nectandra_rbcL_raw)) %>%
    left_join(nectandra_dna, by = "genomic_id") %>%
    left_join(nectandra_specimens, by = "specimen_id") %>%
    assert(is_uniq, genomic_id) %>%
    assert(not_na, genomic_id, specimen_id, taxon, specimen) %>%
    mutate(tip_label = paste(taxon, specimen) %>% 
             str_replace_all(" ", "_")) %>%
    pull(tip_label)
  
  nectandra_rbcL_raw
  
}

#' Rename and align Nectandra rbcL sequences
#' 
#' Renames newly sequenced samples to "Species_Voucher"
#'
#' Also adds sequences of missing taxa from GenBank
#'
#' @param nectandra_rbcL_raw Unaligned rbcL sequences of pteridophytes from Nectandra
#' @param genbank_rbcL_raw Sequences downloaded from GenBank
#'
#' @return DNA alignment of all Nectandra pteridophyte species plus those
#' sequences from GenBank for those that are missing rbcL
#' 
align_nectandra_rbcL <- function (nectandra_rbcL_raw, genbank_rbcL_raw) {
  
  # Combine new sequences with GenBank sequences
  c(nectandra_rbcL_raw, genbank_rbcL_raw) %>%
    # Align sequences with MAFFT, trim ends, and remove any empty cells
    ips::mafft(x = ., exec = "/usr/bin/mafft") %>%
    # Trim ends missing > 50% of sequences
    ips::trimEnds(nrow(.) * 0.5) %>%
    # Delete empty cells
    ips::deleteEmptyCells()
  
}

#' Make a dataframe of minimum interspecific distances
#' 
#' Takes the alignment, calculates raw distances, then
#' outputs the single smallest distance per species
#' (in case of a tie, only a single value is returned)
#' 
#' @param aln DNA alignment (list of class DNAbin)
#' @return tibble
get_min_inter_dist <- function (aln) {
  # Calculate raw pairwise-distances, excluding pairwise missing bases (pairwise.deletion = TRUE)
  ape::dist.dna(aln, model = "raw", pairwise.deletion = TRUE) %>%
    # Convert to tibble
    broom::tidy() %>%
    # Check that there are no self-matches
    assertr::verify(item1 != item2) %>%
    dplyr::mutate_at(dplyr::vars(item1, item2), as.character) %>%
    # Convert to long format (one distance per taxon per row)
    tidyr::pivot_longer(-distance, names_to = "side", values_to = "taxon") %>%
    # Add species (collapses varieties)
    mutate(species = taxastand::sp_name_only(taxon, sep = "_")) %>%
    # Filter to single smallest distance per species (not including multiples for ties)
    dplyr::group_by(species) %>%
    dplyr::arrange(distance) %>%
    dplyr::slice(1) %>%
    dplyr::select(-side) %>%
    dplyr::ungroup() %>%
    # Final checks
    assertr::assert(assertr::not_na, distance, taxon, species) %>%
    assertr::assert(assertr::is_uniq, species) %>%
    assertr::assert(assertr::within_bounds(0,1), distance)
}

#' Bin minimum interspecific distances, with special bin for zeros
#'
#' @param data Tibble of minimum interspecific distances;
#' output of get_min_inter_dist()
#' @param width Bin width (default 0.5%)
#'
#' @return Tibble
#' 
bin_min_inter_dist <- function (data, width = 0.005) {
  
  # Count zero-distances and separate from the rest of the data
  zeros <- data %>% filter(distance == 0) %>% nrow %>%
    tibble(range = 0, n = .)
  
  non_zeroes <-  data %>% filter(distance > 0)
  
  # Bin distances
  ggplot2::cut_width(x = non_zeroes$distance, width = width, boundary = width, closed = "left") %>%
    base::table() %>%
    broom::tidy() %>%
    purrr::set_names(c("range", "n")) %>%
    # Label bins by the right side (upper bound)
    mutate(range = str_split(range, ",") %>% 
             map_chr(2) %>% 
             str_remove("\\)") %>% parse_number) %>% 
    # Add zeros to end
    bind_rows(zeros) %>%
    arrange(range) %>%
    mutate(is_zero = case_when(range == 0 ~ TRUE, TRUE ~ FALSE)) %>%
    mutate(percent = n / sum(n)) %>%
    # Final checks
    assert(not_na, everything()) %>%
    assert(is_uniq, range) %>%
    assert(within_bounds(0, 1), percent)
  
}

#' Get binned minimum interspecific distances from a particular dataset
#'
#' @param full_aln Alignment with species names tagged by their
#' dataset, e.g. "_CR", "_JA"
#' @param dataset_select The dataset to use to calculate
#' minimum interspecific distances
bin_min_inter_dist_by_dataset <- function(full_aln, dataset_select) {
  # Subset the full alignment to only the selected dataset
  full_aln[str_detect(rownames(full_aln), paste0("_", dataset_select, "$")),] %>%
    # Get the minimum interspecific distance for each species
    get_min_inter_dist %>%
    # Bin the distances
    bin_min_inter_dist %>%
    # Add a dataset column
    mutate(dataset = dataset_select)
}

#' Align all rbcL datasets
#'
#' @param nectandra_rbcL rbcL sequences of pteridophytes of Nectandra
#' @param moorea_rbcL rbcL sequences of pteridophytes of Moorea
#' @param japan_rbcL rbcL sequences of pteridophytes of Japan
#' @param japan_rbcL_sexdip rbcL sequences of pteridophytes of Japan,
#' sexual diploids only
#'
#' @return Tibble DNA alignment. Sequences from each dataset will have codes
#' appended to end of name as identifiers: "CR" (Costa Rica, i.e. Nectandra), 
#' "FP" (French Polynesia, i.e. Moorea), 
#' "JA" (Japan), "JAsexdip" (Japan, sexual-diploids only).
#' 
align_all_rbcL <- function(nectandra_rbcL, moorea_rbcL, japan_rbcL, japan_rbcL_sexdip) {
  
  # Add "_CR" to end of name of Nectandra sequences
  nectandra_rbcL <- nectandra_rbcL %>% 
    as.list %>%
    set_names(., paste0(names(.), "_CR"))
  
  # Make sure that all datasets are labeled correctly
  assertthat::assert_that(all(str_detect(names(japan_rbcL), "_JA")))
  assertthat::assert_that(all(str_detect(names(japan_rbcL), "_FP|_CR|JAsexdip", negate = TRUE)))
  assertthat::assert_that(all(str_detect(names(moorea_rbcL), "_FP")))
  assertthat::assert_that(all(str_detect(names(moorea_rbcL), "_JA|_CR|JAsexdip", negate = TRUE)))
  assertthat::assert_that(all(str_detect(names(nectandra_rbcL), "_CR")))
  assertthat::assert_that(all(str_detect(names(nectandra_rbcL), "_JA|_FP|JAsexdip", negate = TRUE)))
  assertthat::assert_that(all(str_detect(names(japan_rbcL_sexdip), "_JAsexdip")))
  assertthat::assert_that(all(str_detect(names(japan_rbcL_sexdip), "_JA$|_FP|_CR", negate = TRUE)))

  # Combine all rbcL sequences
  rbcL_combined <- c(japan_rbcL, moorea_rbcL, nectandra_rbcL, japan_rbcL_sexdip)
  
  # Make global alignment
  rbcL_aln <- mafft(rbcL_combined, exec = "/usr/bin/mafft")
  
  # Trim ends
  trimEnds(rbcL_aln, nrow(rbcL_aln) * 0.5)
  
}

# GenBank accession table ----

#' Make a table of GenBank accession numbers
#'
#' @param new_genbank_accs Newly assigned GenBank accession numbers for rbcL sequences from Nectandra
#' @param nectandra_rbcL_raw Unaligned DNA sequences of newly sequenced specimens from Nectandra;
#' named by genomic accession number ('JNG001', etc)
#' @param DNA_accessions DNA accession numbers and corresponding specimen ID codes
#' @param specimens Specimen data with specimen ID codes
#' @param genbank_rbcL_metadata Metadata of rbcL sequences downloaded from GenBank
#' @param nectandra_rbcL_aligned Aligned Nectandra rbcL sequences used for phylogenetic analysis
#'
#' @return Tibble
#' 
make_genbank_accession_table <- function (new_genbank_accs, nectandra_rbcL_raw, DNA_accessions, specimens, genbank_rbcL_metadata, nectandra_rbcL_aligned) {
  
  # Add scientific name to genbank metadata
  genbank_rbcL_metadata <-
  genbank_rbcL_metadata %>% 
    left_join(
      specimens %>%
        select(taxon, scientific_name) %>%
        unique,
      by = "taxon"
    ) %>%
    assert(not_na, everything())
  
  # Format table of newly sequenced sequences
  tibble(genomic_id = names(nectandra_rbcL_raw)) %>%
    left_join(select(DNA_accessions, genomic_id, specimen_id), by = "genomic_id") %>%
    left_join(select(specimens, specimen_id, taxon, scientific_name, specimen), by = "specimen_id") %>% 
    # Add GenBank accession numbers
    left_join(new_genbank_accs, by = "genomic_id") %>%
    assert(is_uniq, genomic_id, specimen_id) %>%
    assert(not_na, genomic_id, specimen_id, taxon, specimen, genbank_accession) %>%
    # Add "tip_label" in same format as in alignment and tree to check all samples
    # are accounted for
    mutate(tip_label = paste(taxon, specimen) %>% 
             str_replace_all(" ", "_")) %>%
  # Combine with genbank seqs
  bind_rows(genbank_rbcL_metadata) %>%
    # Check that all specimens are identical between GenBank table and alignment
    verify(isTRUE(all.equal(
      sort(tip_label),
      sort(rownames(nectandra_rbcL_aligned))
    ))) %>%
    assert(not_na, taxon, specimen, genbank_accession) %>%
    select(taxon, scientific_name, genomic_id, specimen, genbank_accession) %>%
    arrange(taxon, genomic_id)
  
}

# Etc ----

# Extract a numeric value preceding a keyword in a text string.
get_number <- function(text, keyword) {
  
  str_match(text, glue::glue("([:digit:]+) {keyword}")) %>% magrittr::extract(,2)
  
}

# Plotting ----

#' Define ggplot theme
#' 
#' BW theme with no gridlines, black axis text, main font size 11,
#' axis ticks size 9.
#'
standard_theme2 <- function () {
  ggplot2::theme_bw() + 
    theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour="black"),
      axis.text.y = ggplot2::element_text(colour="black")
    )
}

#' Define ggplot theme
#' 
#' BW theme with no gridlines, black axis text, main font size 11,
#' axis ticks size 9.
#'
standard_theme3 <- function () {
  ggplot2::theme_bw() + 
    theme(
      axis.text.x = ggplot2::element_text(colour="black"),
      axis.text.y = ggplot2::element_text(colour="black")
    )
}

#' Make a richness extrapolation plot
#'
#' @param inext_out Results of iNEXT::iNEXT() on species'
#' occurrences
#'
#' @return ggplot object
make_inext_plot <- function(inext_out) {
  
  # create subset of observed points to add to plot as dots
  observed_data <- inext_out$iNextEst[inext_out$iNextEst$method == "observed",]
  
  # make plot
  ggplot (inext_out$iNextEst, aes(x=t, y=qD)) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), fill = "grey70", alpha=0.2) +
    geom_line(aes(linetype = method), size=0.8) +
    scale_linetype_manual(values = c("dotted", "solid", "solid")) +
    geom_point(data = observed_data, size=2.5) +
    scale_y_continuous("Richness (no. spp.)") +
    scale_x_continuous("Sampling units (days)") +
    standard_theme3() +
    theme(legend.position = "none", 
          legend.title = element_blank())
}

#' Print rbcL tree to pdf
#'
#' @param phy Phylogenetic tree
#' @param ppgi Taxonomy of pteridophytes following Pteridophyte Phylogeny Group I
#' @param outfile Name of file to use to save tree pdf
#'
#' @return Nothing; externally, the tree will be written as a pdf
#'
plot_nectandra_rbcL_tree <- function(phy, ppgi, outfile) {
  
  # Re-root tree on lycophytes and ladderize
  phy <- root_on_lycos(phy, ppgi)
  
  # Reformat tip labels now that tip order has changed
  new_tips <-
    tibble(new_tip = phy$tip.label) %>%
    # Remove "Nitta" and country names from tips
    mutate(
      new_tip = str_remove_all(new_tip, "Nitta_|_Costa_Rica|_Bolivia|_Venezuela") %>%
        str_replace_all("Abrodictyum_rigidum_Dubuisson_HV_1997-3", "Abrodictyum_rigidum_Dubuisson_HV-1997-3"))
  
  phy$tip.label <- new_tips$new_tip
  
  # Reformat node labels
  node_labels_tibble <- tibble(
    old_lab = phy$node.label
  ) %>%
    separate(old_lab, into = c("alrt", "bs"), sep ="\\/", remove = FALSE) %>%
    mutate_at(vars(alrt, bs), ~parse_number(.) %>% round) %>%
    mutate(alrt = case_when(
      alrt < 50 ~ "-",
      alrt == 100 ~ "*",
      TRUE ~ as.character(alrt)
    )) %>%
    mutate(bs = case_when(
      bs < 50 ~ "-",
      bs == 100 ~ "*",
      TRUE ~ as.character(bs)
    )) %>%
    mutate(new_lab = case_when(
      !is.na(alrt) & !is.na(bs) ~ paste(alrt, bs, sep = "/"),
      TRUE ~ old_lab
    ))
  
  phy$node.label <- node_labels_tibble$new_lab
  
  #### Plot using ggtree
  
  # Left side: tree without branch lengths, 
  # include support values and tip labels
  tree1 <- ggtree(phy, branch.length="none") +
    geom_nodelab(hjust = 0, size = 1*.pt) + 
    geom_tiplab() +
    # Need some extra space for long tip names
    theme(plot.margin = margin(t=0,l=0,b=0,r=2.2, unit = "in")) +  
    coord_cartesian(clip = "off")
  
  # Right side: tree with branch lengths, reversed, 
  # no support values or tip labels
  tree2 <- ggtree(phy) + 
    scale_x_reverse() + 
    # Add scale
    geom_treescale(x = -0.1, y = 150, offset = 2) +
    # Make sure margins are set same on both trees so tips line up
    theme(
      plot.margin = margin(t=0,l=0,b=0,r=0, unit = "in"),
      plot.background = element_blank(),
      panel.background = element_blank())
  
  # Combine left and right sides
  final_plot <- tree1 + tree2 + plot_layout(widths = c(2,1))
  
  # Output pdf
  print_height <- 0.16 * length(phy$tip.label)
  
  ggsave(plot = final_plot, filename = outfile, width = 17, height = print_height, units = "in")
  
}

#' Print rbcL tree to pdf
#'
#' @param rbcL_tree Phylogenetic tree
#' @param ppgi Taxonomy of pteridophytes following Pteridophyte Phylogeny Group I
#' @param outfile Name of file to use to save tree pdf
#'
#' @return Nothing; externally, the tree will be written as a pdf
#'
plot_family_rbcL_tree <- function(rbcL_tree, ppgi, outfile) {
  
  # Extract tips into tibble and add taxonomy
  tips <-
    tibble(tip = rbcL_tree$tip.label) %>%
    mutate(
      species = taxastand::sp_name_only(tip, sep = "_"),
      genus = taxastand::genus_name_only(tip, sep = "_")) %>%
    left_join(select(ppgi, genus, family, class))
  
  # Identify lycophyte tips for rooting
  lycos <- tips %>% filter(class == "Lycopodiopsida") %>% pull(tip)
  
  # Root on lycophytes
  rbcL_tree <-
    ape::root(rbcL_tree, outgroup = lycos) %>%
    ape::ladderize(right = FALSE)
  
  # Relabel as families only
  # make tip tabel again now that tree is in different order
  tips <-
    tibble(tip = rbcL_tree$tip.label) %>%
    mutate(
      species = taxastand::sp_name_only(tip, sep = "_"),
      genus = taxastand::genus_name_only(tip, sep = "_")) %>%
    left_join(select(ppgi, genus, family, class))
  
  rbcL_tree$tip.label <- tips$family
  
  # Format bootstrap value: only print BS > 50
  bs <-
    rbcL_tree$node.label %>% 
    parse_number
  
  bs <- case_when(bs > 50 ~ bs) %>% 
    as.character() %>%
    replace_na("")
  
  rbcL_tree$node.label <- bs
  
  # Output pdf
  print_height = 0.114 * length(rbcL_tree$tip.label)
  
  pdf(file = outfile, width = 7, height = print_height)
  plot(rbcL_tree, cex = 0.6, no.margin = TRUE)
  # add node support values
  nodelabels(rbcL_tree$node.label, adj=c(1.2, -0.3), frame="n", cex=0.35, font=1, col="darkgrey")
  # add scale bar
  add.scale.bar(x = 0, y = length(rbcL_tree$tip.label), cex=0.5)
  dev.off()
}

#' Print rbcL tree to pdf
#'
#' @param phy Phylogenetic tree
#' @param outgroup Character vector; tips to use as the outgroup
#' @param outfile Name of file to use to save tree pdf
#' @param nodelab_size Font size to use for node labels
#'
#' @return Nothing; externally, the tree will be written as a pdf
#'
plot_broad_rbcL_tree <- function(phy, outgroup, outfile, nodelab_size = 1) {
  
  # Root tree and ladderize
  phy <-
    ape::root(phy, outgroup = outgroup, resolve.root = TRUE) %>%
    ape::ladderize(right = FALSE)
  
  # Reformat node labels (input tree is from fasttree, so node values are '0.99', etc.)
  node_labels_tibble <- tibble(
    old_lab = phy$node.label
  ) %>%
    # Convert to percent
    mutate(new_lab = parse_number(old_lab) * 100) %>%
    mutate(new_lab = case_when(
      # Don't print for unsupported nodes
      new_lab < 50 ~ "",
      TRUE ~ as.character(new_lab))) %>%
    mutate(new_lab = replace_na(new_lab, ""))
  
  phy$node.label <- node_labels_tibble$new_lab
  
  # Make data frame of data for annotation
  # (just indicating if taxon is a "focus" taxon,
  # i.e, a Nitta collection in the ingroup)
  annotation_data <- tibble(id = phy$tip.label) %>%
    mutate(
      focus_taxon = case_when(
        id %in% outgroup ~ "no",
        str_detect(id, "Nitta") ~ "yes",
        TRUE ~ "no"
      )
    )
  
  #### Plot using ggtree
  
  # Left side: tree without branch lengths, 
  # include support values and tip labels
  tree1 <- ggtree(phy, branch.length="none") %<+% 
    # Merge with data for annotation with special ggtree `%<+%` operator
    annotation_data +
    # First make layer of "labels" at each tip that are just
    # empty boxes highlighted for focus taxon
    geom_tiplab(geom = "label", aes(fill = focus_taxon, color = focus_taxon)) +
    scale_fill_manual(values = c("yes" = "yellow", "no" = "transparent")) +
    scale_color_manual(values = c("yes" = "transparent", "no" = "transparent")) +
    # Add the node labels
    geom_nodelab(hjust = 0, size = nodelab_size*.pt) + 
    # Add the tip labels
    geom_tiplab() +
    # Need some extra space for long tip names
    theme(
      plot.margin = margin(t=0,l=0,b=0,r=3.5, unit = "in"),
      legend.position = "none") +  
    coord_cartesian(clip = "off")
  
  # Right side: tree with branch lengths, reversed, 
  # no support values or tip labels
  tree2 <- ggtree(phy) + 
    scale_x_reverse() + 
    # Add scale
    geom_treescale(x = -0.01, y = length(phy$tip.label) - 1, offset = 2) +
    # Make sure margins are set same on both trees so tips line up
    theme(plot.margin = margin(t=0,l=0,b=0,r=0, unit = "in"))
  
  # Combine left and right sides
  final_plot <- tree1 + tree2 + plot_layout(widths = c(2,1))
  
  # Output pdf
  print_height <- 0.16 * length(phy$tip.label)
  
  ggsave(plot = final_plot, filename = outfile, width = 17, height = print_height, units = "in", limitsize = FALSE)
  
}

# Phylogenetic analysis ----

#' Make a 'broad' DNA alignment including Nectandra rbcL sequences and all
#' rbcL sequences from GenBank for a group of interest
#'
#' @param nectandra_rbcL List of class 'DNAbin'; Nectandra rbcL sequences
#' @param ppgi Dataframe; Pteridophyte phylogeny group I taxonomy
#' @param nectandra_family String; name of family of interest in Nectandra data
#' @param genbank_group String; name of group to download rbcL from GenBank
#' @param start_date Earliest date to include for sequences
#' @param end_date Latest date to include for sequences
#' @param exclude_list Character vector of GenBank accessions to exclude from
#' results (misidentifications)
#'
#' @return List of class 'DNAbin'; Combined Nectandra and GenBank rbcL sequences
#' 
make_broad_alignment <- function(
  nectandra_rbcL, ppgi,
  nectandra_family, genbank_group,
  start_date = "1980/01/01",
  end_date = "2020/03/01",
  exclude_list = "") {
  
  # Extract target family from Nectandra rbcL alignment
  taxa_selected <- 
    tibble(tips = rownames(nectandra_rbcL)) %>%
    mutate(genus = str_split(tips, "_") %>% map_chr(1)) %>%
    left_join(ppgi, by = "genus") %>%
    filter(family %in% nectandra_family) %>%
    pull(tips)
  
  nectandra_rbcL_selected <- nectandra_rbcL[taxa_selected, ] %>% as.list
  
  # Download all rbcL seqs for group of interest from genbank
  
  genbank_rbcL <- gbfetch::fetch_sequences(
    glue::glue('{genbank_group}[ORGN] AND rbcl[Gene] AND 1000:1600[SLEN] NOT accd[Gene] NOT atpB[Gene] NOT spacer AND ("{start_date}"[PDAT]:"{end_date}"[PDAT])'))
  
  genbank_rbcL_metadata <- gbfetch::fetch_metadata(
    glue::glue('{genbank_group}[ORGN] AND rbcl[Gene] AND 1000:1600[SLEN] NOT accd[Gene] NOT atpB[Gene] NOT spacer AND ("{start_date}"[PDAT]:"{end_date}"[PDAT])'))
  
  # Exclude any sequences in exclude list
  genbank_rbcL <- genbank_rbcL[!names(genbank_rbcL) %in% exclude_list]
  genbank_rbcL_metadata <- filter(genbank_rbcL_metadata, !accession %in% exclude_list)
  
  # Rename GB sequences by species + accession
  gb_names <- tibble(accession = names(genbank_rbcL)) %>%
    left_join(select(genbank_rbcL_metadata, accession, species), by = "accession") %>%
    # Truncate extremely long names so they will fit on the tree
    mutate(species_short = str_trunc(species, 40, side = "right", ellipsis = "")) %>%
    mutate(new_name = paste3(species_short, accession) %>% str_replace_all(" ", "_"))
  
  names(genbank_rbcL) <- gb_names$new_name
  
  # Combine and align
  combined_rbcL <- c(nectandra_rbcL_selected, genbank_rbcL)
  
  combined_rbcL_aln <-
    ips::mafft(combined_rbcL, exec = "/usr/bin/mafft", options = "--adjustdirection") %>%
    trimEnds(nrow(.) * 0.5) %>%
    deleteEmptyCells()
  
  # Remove any offending characters: perioids and single quotes
  rownames(combined_rbcL_aln) <- rownames(combined_rbcL_aln) %>% 
    str_remove_all("'") %>% 
    str_remove_all("\\.")

  combined_rbcL_aln

}

#' Identify taxa to use as an outgroup in a tree
#'
#' @param phy Pteridophyte phylogeny. Tip names should be formatted
#' as genus_specific-epithet_accession
#' @param ppgi Pteridophyte Phylogeny Group I taxonomic system
#' @param ... Specification for outgroup taxa used as input to dplyr::filter()
#'
#' @return Character vector
#' 
identify_outgroup <- function(phy, ppgi, ...) {
  
  # Extract tips into tibble and add taxonomy
  tips <-
    tibble(tip = phy$tip.label) %>%
    mutate(
      species = taxastand::sp_name_only(tip, sep = "_"),
      genus = taxastand::genus_name_only(tip, sep = "_")) %>%
    left_join(ppgi, by = "genus") 
  
  # Identify outgroup taxa for rooting
  tips %>% filter(...) %>% pull(tip)
}

#' Root tree on lycophytes
#'
#' Also ladderize
#'
#' @param tree Phylogenetic tree
#' @param ppgi PPGI taxonomic system
#'
#' @return Rooted tree
root_on_lycos <- function(tree, ppgi) {
  
  # Extract tips into tibble and add taxonomy
  tips_for_rooting <-
    tibble(tip = tree$tip.label) %>%
    mutate(
      species = taxastand::sp_name_only(tip, sep = "_"),
      genus = taxastand::genus_name_only(tip, sep = "_")) %>%
    assert(not_na, genus) %>%
    left_join(dplyr::select(ppgi, genus, class), by = "genus") %>%
    assert(not_na, class)
  
  # Identify lycophyte tips for rooting
  lycos <- tips_for_rooting %>% filter(class == "Lycopodiopsida") %>% pull(tip)
  
  # Root on lycophytes
  ape::root(tree, outgroup = lycos) %>%
    ape::ladderize(right = FALSE)
  
}

#' Get formatted node support value for the most
#' recent common ancestor (MRCA) of a set of
#' tips from phylogeny output by IQTREE
#'
#' @param tree Phylogeny made with IQTREE
#' @param tips Character vector of tips
#'
#' @return String
#' 
get_node_support <- function (tree, tips) {
  
  # Make a tibble to match and extract full tip names from input
  tip_df <- tibble(tip = tree$tip.label)
  
  # Match up tip names to input (could just be the species name)
  tips <- filter(tip_df, str_detect(tip, paste(tips, collapse = "|"))) %>%
    # Make sure we don't match more tip labels than input
    verify(nrow(.) == length(tips)) %>%
    pull(tip)
  
  # Make sure the tips specified are in the tree
  assertthat::assert_that(all(tips %in% tree$tip.label))
  
  # Get the MRCA
  mrca_node <- ape::getMRCA(phy = tree, tip = tips)
  
  # Get the node label corresponding to the MRCA
  # https://stackoverflow.com/questions/51696837/r-phylo-object-how-to-connect-node-label-and-node-number
  node_label <- tree$node.label[mrca_node - ape::Ntip(tree)]
  
  # Extract support value and make it pretty
  tibble(node_label = node_label) %>%
    # Since this is a IQTree, we have SH-aLRT and UFboot support
    # values separated by a slash.
    # Split these up
    separate(node_label, into = c("SHaLRT", "UFboot"), sep = "\\/") %>%
    mutate_all(as.numeric) %>%
    # Convert values less than 50% to '<50'
    mutate(SHaLRT = case_when(
      SHaLRT < 50 ~ "<50",
      TRUE ~ as.character(SHaLRT)
    )) %>%
    mutate(UFboot = case_when(
      UFboot < 50 ~ "<50",
      TRUE ~ as.character(UFboot)
    )) %>%
    # Format final printed output
    mutate(
      val = glue::glue("SH-aLRT {SHaLRT}%; UFboot {UFboot}%")
    ) %>%
    mutate(val = as.character(val)) %>%
    pull(val)
}

# MS rendering ----

#' Render a number with SI units
#'
#' @param num Number
#' @param suffix SI suffix to attach to number
#' @param auto_space Logical; should a space automatically be inserted
#' between the number and the suffix (default = TRUE)
#' @param accuracy Decimal points to round to (e.g, 0.1)
#' @param scale A scaling factor: num will be multiplied by scale before formating. 
#' @param prefix Symbol to display before value.
#' @param big.mark Character used between every 3 digits to separate thousands.
#' @param decimal.mark The character to be used to indicate the numeric decimal point
#' @param trim Logical, if FALSE, values are right-justified to a common width (see base::format()).
#' @param ... Other arguments passed on to base::format().
#'
#' @return Character vector
#'
#' @examples
#' si(10000, "m")
si <- function(num, suffix = "", auto_space = TRUE,
               accuracy = NULL, scale = 1, prefix = "",
               big.mark = ",", decimal.mark = ".", trim = TRUE, ...) {
  if(auto_space == TRUE & suffix != "") suffix <- paste("", suffix, collapse = " ")
  if(!is.numeric(num)) num <- readr::parse_number(num)
  scales::comma(num, accuracy = accuracy, scale = scale, prefix = prefix, suffix = suffix,
                big.mark = big.mark, decimal.mark = decimal.mark, trim = trim, ...)
}

# Tweaked version of cairo_ps device function with higher-resolution
# for fallback zones
cairo_ps_high_res <- function(...) {
  cairo_ps(fallback_resolution = 1200, ...)
}

# Dummy function to track arbitary output from rmarkdown::render()
render_tracked <- function (tracked_output, dep1 = NULL, dep2 = NULL, ...) {
  rmarkdown::render(...)
}

#' Convert latex file to docx
#' 
#' Requires pandoc to be installed and on command line
#'
#' @param latex Path to input latex file.
#' @param docx Path to output docx file.
#' @param template Path to template docx file.
#' @param wd Working directory to run conversion. Should be same as
#' directory containing any files needed to render latex to pdf.
#'
#' @return List including STDOUT of pandoc; externally, the
#' docx file will be rendered in `wd`.
#' 
latex2docx <- function (latex, docx, template = NULL, wd = getwd()) {
  
  assertthat::assert_that(assertthat::is.readable(latex))
  
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(docx)))
  
  latex <- fs::path_abs(latex)
  
  docx <- fs::path_abs(docx)
  
  template <- if (!is.null(template)) {
    glue::glue("--reference-doc={fs::path_abs(template)}")
  } else {
    NULL
  }
  
  processx::run(
    command = "pandoc",
    args = c("-s", latex, template, "-o", docx),
    wd = wd
  )
  
}

#' Check the order of figure and table citations in an Rmd file
#'
#' @param rmd_file Path to the Rmd file
#' @param cap_type Type of caption; should match the function name
#' used by captioner to make the citation.
#'
#' @examples
#' check_citation_order("ms/nectandra_pteridos.Rmd", "figure")
check_citation_order <- function(rmd_file, cap_type) {
  
  readr::read_lines(rmd_file) %>% 
    stringr::str_split(" ") %>% 
    unlist %>%
    magrittr::extract(., stringr::str_detect(., cap_type)) %>% 
    str_match(glue::glue("^{cap_type}\\(([^\\()]+)\\)")) %>%
    as_tibble(.name_repair = "universal") %>%
    select(key = 2) %>%
    filter(!is.na(key)) %>%
    unique
  
}

#' Make a path to write out a figure
#' 
#' Will save to selected folder with the file extension,
#' replacing spaces with underscores
#'
#' @param fig_key String; figure name set in captions.R
#' @param ext String; file extension
#' @param folder String; name of folder to write file
#'
#' @return String; file path
#'
#' @examples
#' source("ms/captions.R")
#' fig_path("map", ".tiff")
fig_path <- function (fig_key, ext, folder = "results/ms") {
  
  figure(fig_key) %>% 
    str_replace_all(" ", "_") %>%
    here(folder, .) %>%
    fs::path_ext_set(ext)
  
}

#' Make a path to write out a supplemental figure
#' 
#' Will save to selected folder with the file extension,
#' replacing spaces with underscores
#'
#' @param fig_key String; figure name set in captions.R
#' @param ext String; file extension
#' @param folder String; name of folder to write file
#'
#' @return String; file path
s_fig_path <- function (fig_key, ext, folder = "results/ms") {
  
  s_figure(fig_key) %>% 
    str_replace_all(" ", "_") %>%
    here(folder, .) %>%
    fs::path_ext_set(ext)
  
}

#' Make a path to write out a supplemental table
#' 
#' Will save to selected folder with the file extension,
#' replacing spaces with underscores
#'
#' @param fig_key String; table name set in captions.R
#' @param ext String; file extension
#' @param folder String; name of folder to write file
#'
#' @return String; file path
s_table_path <- function (fig_key, ext, folder = "results/ms") {
  
  s_table(fig_key) %>% 
    str_replace_all(" ", "_") %>%
    here(folder, .) %>%
    fs::path_ext_set(ext)
  
}

# GenBank submission ----

#' Format sequence metadata (for submission to GenBank)
#'
#' @param sample_data Collection data needed for GenBank submission. Links to `dna_data` by `specimen_id`
#' @param dna_data Genomic accession numbers (`genomic_id` and `specimen_id`)
#' @param seqs DNA sequences to be submitted to GenBank in FASTA format. Names must match
#' `genomic_id` in `dna_data`
#' @param adjust_remainder: Named numeric vector with names ('0', '1', and '2'): the frameshift to use
#' when the initial gap remainder is 0, 1, or 2. For example, `c("0" = 2, "1" = 1, "2" = 3)`
#' @param manual_readframe Named character vector; manually-specified readframe. Name must match `genomic_id`
#' @param notes Named character vector; notes for specific samples. Name must match `genomic_id` (not currently implemented)
#'
#' @return Dataframe
#'
format_data_for_genbank <- function(sample_data, dna_data, seqs, adjust_remainder = c("0" = 2, "1" = 1, "2" = 3), manual_readframe = NULL, notes = NULL) {
  
  sample_data <-
    sample_data %>%
    assert(not_na, specimen, specimen_id, genus, specific_epithet, herbaria, date_collected, country, locality) %>%
    mutate(date_collected = format(date_collected, "%d-%b-%Y")) %>%
    mutate(
      specimen = str_replace_all(specimen, "Nitta", "J.H. Nitta"),
      certainty = ifelse(is.na(certainty), certainty, paste0(certainty, "."))
    ) %>%
    transmute(
      specimen_id = specimen_id,
      voucher = glue::glue("{specimen} ({herbaria})"),
      country = glue::glue("{country}: {locality}"),
      collection_date = date_collected,
      organism = jntools::paste3(genus, certainty, specific_epithet, infraspecific_rank, infraspecific_name, sep = " "),
      # remove accents, etc from authority, which are not allowed by genbank
      authority = stringi::stri_trans_general(author, "Latin-ASCII"),
    ) %>%
    left_join(dna_data, by = "specimen_id") %>%
    verify(all(names(seqs) %in% genomic_id)) %>%
    filter(genomic_id %in% names(seqs))
  
  # Check that the adjust_remainder argument is properly formatted
  assertthat::assert_that(all(names(adjust_remainder) == c("0", "1", "2")))
  assertthat::assert_that(is.numeric(adjust_remainder))
  
  # Align sequences
  aligned_seqs <- ips::mafft(seqs, exec = system("which mafft", intern = TRUE))
  
  # Get sequence lengths
  seq_lengths <-
    seqs %>%
    map_dbl(length) %>% 
    tibble(genomic_id = names(.), seq_len = .)
  
  # Calculate codon start position based on gap remainder 
  # (divide number of initial gaps by 3, take remainder)
  codon_start_positions <-
    # Convert alignment to a dataframe with
    # rownames as seq ID and a single column with the sequence
    aligned_seqs %>% 
    as.character() %>%
    t %>%
    as.data.frame() %>%
    map(~paste(., collapse = "")) %>%
    # Count initial gaps
    map(~str_match(., "^(-*)") %>% magrittr::extract(,2)) %>%
    map_dbl(~str_count(., "-")) %>%
    tibble(genomic_id = names(.), num_gaps = .) %>%
    # Get remainder after dividing by 3
    mutate(remainder = num_gaps %% 3) %>%
    # Add codon start position
    mutate(codon_start = case_when(
      remainder == 0 ~ adjust_remainder[["0"]],
      remainder == 1 ~ adjust_remainder[["1"]],
      remainder == 2 ~ adjust_remainder[["2"]]
    )) %>%
    select(genomic_id, codon_start)
  
  # Adjust readframe manually if needed
  if (!(is.null(manual_readframe))) {
    
    # Convert input to tibble
    manual_readframe <-
      tibble(
        materialSampleID = names(manual_readframe),
        codon_start = manual_readframe
      )
    
    # Replace with manual codon start position
    codon_start_positions <-
      codon_start_positions %>%
      anti_join(manual_readframe, by = "materialSampleID") %>%
      bind_rows(manual_readframe)
  }
  
  # Convert notes to tibble
  # notes <- tibble(
  #   materialSampleID = names(notes),
  #   note = notes
  # )
  
  # Combine metadata
  combined_sample_data <-
    sample_data %>%
    left_join(codon_start_positions, by = "genomic_id") %>%
    left_join(seq_lengths, by = "genomic_id") %>%
    # Make sure there are no duplications
    assert(is_uniq, genomic_id) %>%
    # Only `authority` and `notes` are allowed to be missing 
    # (authority, in case of taxa only identified to genus)
    assert(not_na, genomic_id, voucher, country, organism, collection_date)
  
  # Attach metadata to sequence names (same as `materialSampleID`)
  tibble(genomic_id = rownames(aligned_seqs)) %>%
    # Add collection data for all samples
    left_join(combined_sample_data, by = "genomic_id") %>%
    # Make sure that worked as planned
    assert(is_uniq, genomic_id) %>%
    assert(not_na, genomic_id, voucher, country, organism)
  
}

#' Format FASTA sequences for submission to GenBank via tbl2asn
#'
#' The metadata will be inserted into the sequence header as described at
#' https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/ and
#' https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html
#'
#' @param seqs_data Sequence metadata
#' @param seqs Sequences in FASTA format
#'
#' @return List of class 'DNAbin'
#'
format_fasta_for_tbl2asn <- function (seqs_data, seqs) {
  
  # Define helper functions to add genbank modifiers to a variable
  # - for a single string
  add_genbank_modifiers_single <- function (x, mod_name) {
    # Return NA if the input is NA
    if(is.na(x)) return(NA)
    glue::glue("[{mod_name}={x}]")
  }
  # - vectorized
  add_genbank_modifiers <- function (x, mod_name) {
    map_chr(x, ~add_genbank_modifiers_single(., mod_name))
  }
  
  # Format names for tbl2asn
  tbl2asn_labels_df <-
    seqs_data %>%
    transmute(
      genomic_id = genomic_id,
      voucher = add_genbank_modifiers(voucher, "specimen-voucher"),
      country = add_genbank_modifiers(country, "country"),
      collection_date = add_genbank_modifiers(collection_date, "collection-date"),
      organism = add_genbank_modifiers(organism, "organism"),
      authority = add_genbank_modifiers(authority, "authority"),
      organelle = "[location=chloroplast]", # this will result in 'organelle="plastid:chloroplast" in the flat-file'
      gcode = "[gcode=11]" # this will result in 'transl_table=11' in the flat-file
    ) %>%
    # Use jntools::paste3() to avoid NAs in the output 
    # e.g., for some rows collection_date might be NA, but this way we won't get an "NA" in the label
    mutate(label = jntools::paste3(
      genomic_id, voucher, country, collection_date, organism, authority, organelle, gcode)) %>%
    select(genomic_id, label)
  
  tbl2asn_labels <-
    tibble(genomic_id = names(seqs)) %>%
    left_join(tbl2asn_labels_df, by = "genomic_id") %>%
    assert(not_na, label, genomic_id) %>%
    assert(is_uniq, label, genomic_id) %>%
    magrittr::extract2("label")
  
  # Rename sequences with formatted labels
  names(seqs) <- tbl2asn_labels
  
  seqs
  
}

#' Make an entry in a feature table for tbl2asn
#'
#' Assumes a single gene/CDS for the entire length of the sequence,
#' with no gaps or introns
#' 
#' @param name Sequence name
#' @param seq_len Sequence length
#' @param codon_start Codon start position
#' @param product Name of gene product
#' (default = "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit")
#' @param gene Name of gene
#' (default = "rbcL")
#' @param transl_table RNA translation table
#' (default = 11, bacterial)
#'
#' @return Text formatted as an entry in a tbl2asn feature table
#' 
make_feature <- function (
  name, seq_len, codon_start,
  product = "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit",
  gene = "rbcL",
  transl_table = 11) {
  
  start = 1
  
  glue::glue(">Feature lcl|{name}
<{start}\t>{seq_len}\tgene
\t\t\tgene\t{gene}
<{start}\t>{seq_len}\tCDS
\t\t\tproduct\t{product}
\t\t\tcodon_start\t{codon_start}
\t\t\ttransl_table\t{transl_table}
\t\t\tprotein_id\tlcl|{name}_1")
  
}

#' Run tbl2asn
#' 
#' Runs tbl2asn for a set of rbcL sequences.
#' 
#' See https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/#src
#' 
#' It is tempting to use a CSV of metadata for each sample, but this fails for 
#' certain columns ('gcode', 'location'), so encode this information in the
#' fasta headers instead.
#'
#' @param genbank_template Path to GenBank template generated by filling out
#' form and downloading from: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
#' @param genbank_features Character vector: sequence features.
#' (see https://www.ncbi.nlm.nih.gov/Sequin/table.html)
#' @param seqs List of class 'DNAbin'; DNA sequences (multi-fasta). Names of each entry
#' must contain metadata in brackets as described in tbl2asn docs.
#' @param submission_name String; name to use for submission
#' @param results_dir Path to directory to write output
#' @param ... Extra arguments; not used by this function, but meant for tracking by drake
#'
#' @return A hash of the genbank flat-file. Externally, the following files will
#' be written to `results_dir` (* is `submission_name`):
#' *.sqn: GenBank sqn file
#' *.gbf: GenBank flatfile
#' *.val: Validation report
#' *_discrepancy_report.txt: Discrepancy report
#' 
tbl2asn <- function (genbank_template_file, genbank_features, seqs, submission_name, results_dir, ...) {
  
  # Create temporary working directory
  fs::dir_create(tempdir(), "tbl2asn")
  temp_dir <- fs::path(tempdir(), "tbl2asn")  
  
  # Write out sequences
  ape::write.FASTA(seqs, fs::path(temp_dir, glue::glue("{submission_name}.fsa")))
  
  # Copy submission template
  fs::file_copy(genbank_template_file, fs::path(temp_dir, glue::glue("{submission_name}.sbt")))
  
  # Write out feature table
  readr::write_lines(genbank_features, fs::path(temp_dir, glue::glue("{submission_name}.tbl")))
  
  # Setup arguments for tbl2asn
  tbl2asn_args = c(
    "-a", "s", # specify input as multi-fasta file 
    "-p", ".", # all input files (.fsa, .sbt, and .src) are in working directory
    "-V", "vb", # run verification and output flatfile
    "-t", fs::path(temp_dir, glue::glue("{submission_name}.sbt")),
    "-Z", glue::glue("{submission_name}_discrepancy_report.txt")
  )
  
  # Run tbl2asn
  processx::run("tbl2asn", tbl2asn_args, wd = temp_dir, echo = TRUE)
  
  # Report errors and warnings
  read_lines(fs::path(temp_dir, "errorsummary.val")) %>%
    stringr::str_trim() %>%
    c("tbl2asn messages:", .) %>%
    print()
  
  # Copy output to results directory
  fs::file_copy(
    fs::path(temp_dir, glue::glue("{submission_name}_discrepancy_report.txt")), 
    fs::path(results_dir, glue::glue("{submission_name}_discrepancy_report.txt")),
    overwrite = TRUE)
  
  fs::file_copy(
    fs::path(temp_dir, glue::glue("{submission_name}.gbf")), 
    fs::path(results_dir, glue::glue("{submission_name}.gbf")),
    overwrite = TRUE)
  
  fs::file_copy(
    fs::path(temp_dir, glue::glue("{submission_name}.sqn")), 
    fs::path(results_dir, glue::glue("{submission_name}.sqn")),
    overwrite = TRUE)
  
  fs::file_copy(
    fs::path(temp_dir, glue::glue("{submission_name}.val")), 
    fs::path(results_dir, glue::glue("{submission_name}.val")),
    overwrite = TRUE)
  
  # Cleanup
  hash <- digest::digest(fs::path(temp_dir, glue::glue("{submission_name}.gbf")))
  
  fs::dir_delete(temp_dir)
  
  # Return hash of gbf file
  hash
  
}

tracked_file <- function (path, ...) {
  path
}

#' Extract a translated amino acid sequence from a single entry
#' in a genbank flatfile
#'
#' @param gb_entry String (character vector of length 1);
#' single entry from a genbank flatfile
#'
#' @return amino acid sequence as a character vector, named
#' for the species + voucher
#' 
extract_translation <- function (gb_entry) {
  
  translation <- 
    gb_entry %>%
    paste(sep = "") %>%
    str_remove_all("\n") %>%
    str_remove_all('\"') %>%
    str_match('translation=(.+)ORIGIN') %>%
    magrittr::extract(,2) %>%
    str_remove_all(" ")
  
  voucher <-
    gb_entry %>%
    paste(sep = "") %>%
    str_remove_all("\n") %>%
    str_remove_all('\"') %>%
    str_match('LOCUS +([^ ]+) +') %>%
    magrittr::extract(,2)
  
  species <-
    gb_entry %>%
    paste(sep = "") %>%
    str_remove_all("\n") %>%
    str_remove_all('\"') %>%
    str_match('ORGANISM +(.+) Unclassified') %>%
    magrittr::extract(,2) %>%
    stringr::str_trim(side = "both")
  
  names(translation) <- paste(species, voucher)
  
  translation
}

#' Parse a genbank file and extract all the amino acid sequences
#'
#' @param gbff_path Path to genbank flat file
#'
#' @return List of class "AAbin"
#' 
parse_aa_from_flatfile <- function (gbff_path) {
  read_file(gbff_path) %>%
    # '\\' is delimiter between entries
    str_split("\\/\\/") %>%
    unlist %>%
    # Drop the last item, as it is just an empty line (after the last '\\')
    magrittr::extract(-length(.)) %>%
    # Extract AA sequences from each entry
    map(extract_translation) %>%
    # Name them as the species + voucher
    set_names(map_chr(., names)) %>%
    # Convert to ape format
    ape::as.AAbin()
}
