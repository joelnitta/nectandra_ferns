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
      parse_names_batch(.$scientific_name) %>% 
        select(scientific_name = b, taxon = c),
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
      parse_names_batch(.$scientific_name) %>% 
        select(scientific_name = b, taxon = c),
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
    filter(locality == "Nectandra Cloud Forest Preserve") %>%
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
    filter(locality == "Nectandra Cloud Forest Preserve") %>%
    assert(not_na, date_collected) %>%
    # Convert to presence-absence matrix of species x dates
    select(taxon, date_collected) %>%
    unique %>%
    mutate(abun = 1) %>%
    pivot_wider(names_from = date_collected, values_from = abun, 
                values_fill = list(abun = 0)) %>%
    column_to_rownames("taxon") %>% 
    as.matrix %>%
    iNEXT(datatype = "incidence_raw", endpoint = endpoint)

  }

# Barcode analysis ----

#' Align Nectandra rbcL sequences
#'
#' Also checks for missing taxa and adds them from GenBank
#'
#' @param nectandra_rbcL_raw Unaligned rbcL sequences of pteridophytes from Nectandra
#' @param nectandra_dna DNA accession numbers linking accession to specimen ID
#' @param nectandra_specimens Specimen collection data including species names
#'
#' @return DNA alignment of all Nectandra pteridophyte species plus those
#' sequences from GenBank for those that are missing rbcL
#' 
align_rbcL <- function (nectandra_rbcL_raw, nectandra_dna, nectandra_specimens) {
  
  # Check for missing taxa from Nectandra list.
  nectandra_rbcL_missing_taxa <-
    tibble(genomic_id = names(nectandra_rbcL_raw)) %>%
    left_join(nectandra_dna) %>%
    left_join(select(nectandra_specimens, specimen_id, taxon)) %>%
    assert(not_na, everything()) %>%
    anti_join(nectandra_specimens, ., by = "taxon") %>%
    pull(taxon) %>%
    unique %>% sort
  
  # Construct query and fetch GenBank sequences for missing taxa
  query <- glue::glue("({paste(nectandra_rbcL_missing_taxa, collapse = '[ORIGIN] OR ')}[ORIGIN]) AND rbcl")
  
  gb_seqs <- gbfetch::fetch_sequences(query)
  
  gb_seqs_metadata <- gbfetch::fetch_metadata(query)
  
  # Select the appropriate sequences to use and rename them.
  # KM008147 = Pteris altissima
  # AY175795 = Trichomanes polypodioides
  gb_seqs <- gb_seqs[names(gb_seqs) %in% c("KM008147", "AY175795")]
  
  names(gb_seqs) <-
    tibble(accession = names(gb_seqs)) %>%
    left_join(gb_seqs_metadata) %>%
    mutate(tip_label = paste(species, accession) %>% str_replace_all(" ", "_")) %>%
    pull(tip_label)
  
  # Reassign names to be taxon plus DNA accession number
  names(nectandra_rbcL_raw) <-
    tibble(genomic_id = names(nectandra_rbcL_raw)) %>%
    left_join(nectandra_dna, by = "genomic_id") %>%
    left_join(nectandra_specimens, by = "specimen_id") %>%
    assert(is_uniq, genomic_id) %>%
    assert(not_na, genomic_id, specimen_id, taxon) %>%
    mutate(tip_label = paste(taxon, genomic_id) %>% str_replace_all(" ", "_")) %>%
    pull(tip_label)
  
  # Combine new sequences with GenBank sequences
  nectandra_rbcL_raw <- c(nectandra_rbcL_raw, gb_seqs)
  
  # Align sequences with MAFFT, trim ends, and remove any empty cells
  ips::mafft(nectandra_rbcL_raw, exec = "/usr/bin/mafft") %>%
    trimEnds(nrow(.) * 0.5) %>%
    deleteEmptyCells()
  
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
  ape::dist.dna(aln, model = "raw", pairwise.deletion = TRUE) %>%
    broom::tidy() %>%
    dplyr::mutate_at(dplyr::vars(item1, item2), as.character) %>%
    tidyr::gather(side, taxon, -distance) %>%
    mutate(species = sp_name_only(taxon, sep = "_")) %>%
    dplyr::group_by(species) %>%
    dplyr::arrange(distance) %>%
    dplyr::slice(1) %>%
    dplyr::select(-side) %>%
    dplyr::ungroup()
}

#' Bin minimum interspecific distances, with special bin for zeros
#'
#' @param data Tibble of minimum interspecific distances;
#' output of get_min_inter_dist()
#' @param width Bin width
#'
#' @return Tibble
#' 
bin_min_inter_dist <- function (data, width = 0.005) {
  
  zeros <- data %>% filter(distance == 0) %>% nrow %>%
    tibble(range = 0, n = .)
  
  non_zeroes <-  data %>% filter(distance > 0)
  
  cut_width(non_zeroes$distance, width = width, boundary = width, closed = "left") %>%
    table %>%
    tidy %>%
    set_names(c("range", "n")) %>%
    mutate(range = str_split(range, ",") %>% 
             map_chr(2) %>% 
             str_remove("\\)") %>% str_remove("]")) %>% 
    mutate(range = as.numeric(range)) %>%
    bind_rows(zeros) %>%
    mutate(is_zero = case_when(range == 0 ~ TRUE, TRUE ~ FALSE)) %>%
    mutate(percent = n / sum(n))
  
}

#' Get binned minimum interspecific distances from a particular dataset
#'
#' @param full_aln Alignment with species names tagged by their
#' dataset, e.g. "_CR", "_JA"
#' @param dataset_select The dataset to use to calculate
#' minimum interspecific distances
bin_min_inter_dist_by_dataset <- function(full_aln, dataset_select) {
  full_aln[str_detect(rownames(full_aln), paste0("_", dataset_select, "$")),] %>%
    get_min_inter_dist %>%
    bin_min_inter_dist %>%
    mutate(dataset = dataset_select)
}

#' Calculate minimum interspecific distances across three
#' rbcL datasets
#'
#' @param nectandra_rbcL rbcL sequences of pteridophytes of Nectandra
#' @param moorea_rbcL rbcL sequences of pteridophytes of Moorea
#' @param japan_rbcL rbcL sequences of pteridophytes of Japan
#' @param japan_rbcL_sexdip rbcL sequences of pteridophytes of Japan,
#' sexual diploids only
#'
#' @return Tibble
#' 
analyze_min_dist <- function(nectandra_rbcL, moorea_rbcL, japan_rbcL, japan_rbcL_sexdip) {
  
  # Add "_CR" to end of name of Nectandra sequences
  nectandra_rbcL <- nectandra_rbcL %>% 
    as.list %>%
    set_names(., paste0(names(.), "_CR"))
  
  # Combine all rbcL sequences
  rbcL_combined <- c(japan_rbcL, moorea_rbcL, nectandra_rbcL, japan_rbcL_sexdip)
  
  # Make global alignment
  rbcL_aln <- mafft(rbcL_combined, exec = "/usr/bin/mafft")
  
  # Trim ends
  rbcL_aln <- trimEnds(rbcL_aln, nrow(rbcL_aln) * 0.5)
  
  # Calculate minimum interspecific distances for each
  # dataset separately and combine these into a single dataframe
  map_df(c("CR", "FP", "JA", "JAsexdip"), 
         ~bin_min_inter_dist_by_dataset(full_aln = rbcL_aln, dataset_select = .))
}

# GenBank accession table ----

#' Make a table of GenBank accession numbers
#'
#' @param nectandra_rbcL rbcL alignment used for phylogenetic analysis
#' @param DNA_accessions DNA accession numbers and corresponding specimen ID codes
#' @param specimens Specimen data
#'
#' @return Tibble
#' 
make_genbank_accession_table <- function (nectandra_rbcL, DNA_accessions, specimens) {
  
  tibble(
    tip = rownames(nectandra_rbcL)
  ) %>%
    mutate(genomic_id = str_match(tip, "_([:upper:]+.+)$") %>% magrittr::extract(,2)) %>%
    left_join(select(DNA_accessions, genomic_id, specimen_id), by = "genomic_id") %>%
    # Manually set species ID for Pteris_altissima_KM008147 (Nitta 863)
    mutate(specimen_id = case_when(
      tip == "Pteris_altissima_KM008147" ~ 892,
      TRUE ~ specimen_id
    )) %>%
    left_join(select(specimens, specimen_id, taxon, scientific_name, specimen), by = "specimen_id") %>% 
    mutate(
      genbank_accession = case_when(
        str_detect(genomic_id, "JNG") ~ "TBD",
        TRUE ~ genomic_id
      ),
      genomic_id = case_when(
        str_detect(genomic_id, "JNG") ~ genomic_id,
        TRUE ~ NA_character_
      ),
      # Manually add data for AY175795
      taxon = case_when(
        genbank_accession == "AY175795" ~ "Trichomanes polypodiodes",
        TRUE ~ taxon
      ),
      scientific_name = case_when(
        genbank_accession == "AY175795" ~ "Trichomanes polypodiodes L.",
        TRUE ~ scientific_name
      ),
      specimen = case_when(
        genbank_accession == "AY175795" ~ "M. Kessler 8808 (Bolivia)",
        TRUE ~ specimen
      )
    ) %>%
    select(taxon, scientific_name, genomic_id, specimen, genbank_accession) %>%
    assert(not_na, taxon, specimen, genbank_accession)
  
}

# Etc ----

#' Parse species names in batch
#' 
#' Runs much faster for a large number of names.
#' 
#' Requires gnparser to be installed and on $PATH
#'
#' @param names Character vector of species names. 
#' May include author, variety, etc. All names must be
#' unique, with no NAs.
#' @param check Logical; should a check be made that the
#' results original name match the names of the input?
#'
#' @return Tibble
#'
#' @examples
#' parse_names_batch("Amaurorhinus bewichianus (Wollaston,1860) (s.str.)")
parse_names_batch <- function (names, check = TRUE) {
  
  assertthat::assert_that(is.character(names))
  assertthat::assert_that(all(assertr::not_na(names)))
  assertthat::assert_that(all(assertr::is_uniq(names)))
  
  temp_file <- fs::file_temp() %>% fs::path_ext_set("txt")
  
  temp_dir <- fs::path_dir(temp_file)
  
  temp_txt <- fs::path_file(temp_file)
  
  readr::write_lines(names, temp_file)
  
  args = c(
    "-f",
    "simple",
    "-j",
    "20",
    temp_txt
  )
  
  results <- 
    processx::run(
      command = "gnparser", args, wd = temp_dir) %>%
    magrittr::extract("stdout") %>%
    unlist() %>%
    read_lines() %>%
    stringr::str_split("\\|") %>% 
    # Need to figure out what each column actually means
    purrr::map(~purrr::set_names(., letters[1:7])) %>%
    do.call(bind_rows, .) %>%
    select(b:g)
  
  if (isTRUE(check)) {
    # b contains the original name. Sometimes these can get out of order.
    # Make sure they are in the same order as the input.
    results <-
      arrange(results, match(b,names)) %>%
      assertr::verify(all(.$b == names))
  }
  
  results
  
}


#' Extract only the species name from a longer name.
#' 
#' It is assumed that the first two parts of the name are genus then
#' specific epithet. No checking is done for this.
#'
#' @param taxon_name Taxon name, e.g. "Crepidomanes minutum var minutum".
#' @param sep Character separating parts of the name.
#'
#' @return The first two parts of the name separated by space.
#' @examples
#' sp_name_only("Crepidomanes minutum var minutum")
sp_name_only <- function (taxon_name, sep = " ") {
  
  assertthat::assert_that(is.character(taxon_name))
  
  str_split(taxon_name, sep) %>% 
    map_chr(., ~magrittr::extract(., 1:2) %>% jntools::paste3(collapse = sep))
}

#' Extract only the genus name from a longer name.
#' 
#' It is assumed that the first part of the name is the genus.
#' No checking is done for this.
#'
#' @param taxon_name Taxon name, e.g. "Crepidomanes minutum var minutum".
#' @param sep Character separating parts of the name.
#'
#' @return The first two parts of the name separated by space.
#' @examples
#' genus_name_only("Crepidomanes minutum var minutum")
genus_name_only <- function (taxon_name, sep = " ") {
  
  assertthat::assert_that(is.character(taxon_name))
  
  str_split(taxon_name, sep) %>% 
    map_chr(., ~magrittr::extract(., 1))
}

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
    geom_line(aes(linetype=method), size=0.8) +
    scale_linetype_manual(values=c("dotted", "solid", "solid")) +
    geom_point(data=observed_data, size=2.5) +
    scale_y_continuous("Richness (no. spp.)") +
    scale_x_continuous("Sampling units (days)") +
    standard_theme3() +
    theme(legend.position="none", 
          legend.title=element_blank())
}

#' Print rbcL tree to pdf
#'
#' @param phy Phylogenetic tree
#' @param ppgi Taxonomy of pteridophytes following Pteridophyte Phylogeny Group I
#' @param specimens Specimen data
#' @param outfile Name of file to use to save tree pdf
#'
#' @return Nothing; externally, the tree will be written as a pdf
#'
plot_rbcL_tree <- function(phy, ppgi, specimens, dna_acc, outfile) {
  
  # Extract tips into tibble and add taxonomy
  tips <-
    tibble(tip = phy$tip.label) %>%
    mutate(
      species = sp_name_only(tip, sep = "_"),
      genus = genus_name_only(tip, sep = "_")) %>%
    left_join(dplyr::select(ppgi, genus, family, class), by = "genus") %>%
    mutate(genomicID = str_match(tip, "_([:upper:]+.+)$") %>% magrittr::extract(,2))
  
  # Identify lycophyte tips for rooting
  lycos <- tips %>% filter(class == "Lycopodiopsida") %>% pull(tip)
  
  # Root on lycophytes
  phy <-
    ape::root(phy, outgroup = lycos) %>%
    ape::ladderize(right = FALSE)
  
  # Reformat tip labels now that tip order has changed
  new_tips <-
    tibble(tip = phy$tip.label) %>%
    mutate(genomic_id = str_match(tip, "_([:upper:]+.+)$") %>% magrittr::extract(,2)) %>%
    left_join(dplyr::select(dna_acc, genomic_id, specimen_id), by = "genomic_id") %>%
    # Manually set species ID for Pteris_altissima_KM008147 (Nitta 863)
    mutate(specimen_id = case_when(
      tip == "Pteris_altissima_KM008147" ~ 892,
      TRUE ~ specimen_id
    )) %>%
    left_join(dplyr::select(specimens, specimen_id, taxon, specimen), by = "specimen_id") %>%
    mutate(
      new_tip = case_when(
        !is.na(specimen_id) ~ paste(taxon, specimen),
        TRUE ~ tip)
    ) %>%
    mutate(new_tip = new_tip %>%
             str_remove_all("Nitta ") %>% 
             str_replace_all("_", " ")
    )
  
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
    theme(plot.margin = margin(t=0,l=0,b=0,r=2.5, unit = "in")) +  
    coord_cartesian(clip = "off")
  
  # Right side: tree with branch lengths, reversed, 
  # no support values or tip labels
  tree2 <- ggtree(phy) + 
    scale_x_reverse() + 
    # Add scale
    geom_treescale(x = -0.1, y = 150, offset = 5) +
    # Make sure margins are set same on both trees so tips line up
    theme(plot.margin = margin(t=0,l=0,b=0,r=0, unit = "in"))
  
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
plot_rbcL_tree_families <- function(rbcL_tree, ppgi, outfile) {
  
  # Extract tips into tibble and add taxonomy
  tips <-
    tibble(tip = rbcL_tree$tip.label) %>%
    mutate(
      species = sp_name_only(tip, sep = "_"),
      genus = genus_name_only(tip, sep = "_")) %>%
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
      species = sp_name_only(tip, sep = "_"),
      genus = genus_name_only(tip, sep = "_")) %>%
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

# Phylogenetic analysis ----

#' Make a 'broad' tree including Nectandra rbcL sequences and all
#' rbcL sequences from GenBank for a group of interest
#'
#' @param nectandra_rbcL List of class 'DNAbin'; Nectandra rbcL sequences
#' @param ppgi Dataframe; Pteridophyte phylogeny group I taxonomy
#' @param nectandra_family String; name of family of interest in Nectandra data
#' @param genbank_group String; name of group to download rbcL from GenBank
#' @param tree_path Path to write out phylogenetic tree.
#'
#' @return List; phylogenetic tree in ape format. Externally, the tree
#' will be written in NEXUS format to tree_path
#' 
make_broad_tree <- function(
  nectandra_rbcL, ppgi,
  nectandra_family, genbank_group,
  tree_path) {
  
  # Extract target family from Nectandra rbcL alignment
  taxa_selected <- 
    tibble(tips = rownames(nectandra_rbcL)) %>%
    mutate(genus = str_split(tips, "_") %>% map_chr(1)) %>%
    left_join(ppgi, by = "genus") %>%
    filter(family == nectandra_family) %>%
    pull(tips)
  
  nectandra_rbcL_selected <- nectandra_rbcL[taxa_selected, ] %>% as.list
  
  # Download all rbcL seqs for group of interest from genbank
  
  genbank_rbcL <- gbfetch::fetch_sequences(
    glue::glue("{genbank_group}[ORGN] AND rbcl[Gene] AND 1000:1600[SLEN] NOT accd[Gene] NOT atpB[Gene] NOT spacer"))
  
  genbank_rbcL_metadata <- gbfetch::fetch_metadata(
    glue::glue("{genbank_group}[ORGN] AND rbcl[Gene] AND 1000:1600[SLEN] NOT accd[Gene] NOT atpB[Gene] NOT spacer"))
  
  # Rename GB sequences by species + accession
  gb_names <- tibble(accession = names(genbank_rbcL)) %>%
    left_join(select(genbank_rbcL_metadata, accession, species), by = "accession") %>%
    mutate(new_name = paste3(species, accession) %>% str_replace_all(" ", "_"))
  
  names(genbank_rbcL) <- gb_names$new_name
  
  # Combine and align
  combined_rbcL <- c(nectandra_rbcL_selected, genbank_rbcL)
  
  combined_rbcL_aln <-
    ips::mafft(combined_rbcL, exec = "/usr/local/bin/mafft", options = "--adjustdirection") %>%
    trimEnds(nrow(.) * 0.5) %>%
    deleteEmptyCells()
  
  # Remove any offending characters: perioids and single quotes
  rownames(combined_rbcL_aln) <- rownames(combined_rbcL_aln) %>% str_remove_all("'")
  
  rownames(combined_rbcL_aln) <- rownames(combined_rbcL_aln) %>% str_remove_all("\\.")
  
  # Make tree with fasttree
  combined_rbcL_tree <- fasttree(combined_rbcL_aln)
  
  ape::write.tree(combined_rbcL_tree, tree_path)
  
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
render_tracked <- function (tracked_output, ...) {
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
