# Data wrangling ----

#' Download Nitta et al 2017 Ecol Mono data zip file and 
#' extract needed data files
#'
#' @param dl_path Name of file to download the zip file to.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return Unzipped data files:
#' - rbcL_clean_sporos.fasta: rbcL sequences of sporophytes from Moorea
#'
download_and_unzip_nitta_2017 <- function (dl_path, unzip_path, ...) {
  
  # Make sure the target directory exists
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(dl_path)))
  
  # Set url
  url <- "https://datadryad.org/bitstream/handle/10255/dryad.132050/data_and_scripts.zip?sequence=1"
  
  # Download zip file
  download.file(url, dl_path)
  
  # Unzip only needed data files to data/nitta_2017/
  unzip(dl_path, "data_and_scripts/shared_data/rbcL_clean_sporos.fasta", exdir = unzip_path, junkpaths = TRUE)
  
}

tidy_taxonomy <- function (taxonomy_data) {
  
  taxonomy_data %>%
    select(-genus_species) %>%
    rename(
      specific_epipthet = species,
      infraspecific_epipthet = infrasp_name) %>%
    mutate(taxon = jntools::paste3(genus, specific_epipthet, infraspecific_epipthet)) %>%
    mutate(sci_name = jntools::paste3(genus, specific_epipthet, author, 
                             infrasp_rank, infraspecific_epipthet, var_author)) %>%
    mutate(sci_name = str_trim(sci_name, "both")) %>%
    select(taxon, sci_name)
  
}

tidy_specimens <- function (specimen_data, ppgi, taxonomy) {
  
  specimen_data %>%
    filter(country == "Costa Rica") %>%
    as_tibble %>%
    filter(locality %in% c(
      "Nectandra Cloud Forest Preserve",
      "Finca Ocotea"
    )) %>%
    filter(is_gametophyte == 0) %>%
    filter(!is.na(species)) %>%
    rename(specific_epipthet = species, infraspecific_epipthet = infraspecific_name) %>%
    mutate(
      species = paste(genus, specific_epipthet),
      taxon = paste3(genus, specific_epipthet, infraspecific_epipthet)) %>%
    left_join(taxonomy) %>%
    left_join(select(ppgi, genus, family, class)) %>%
    mutate(sci_name = case_when(
      is.na(sci_name) ~ taxon,
      TRUE ~ sci_name
    )) %>%
    select(class, family, genus, species, taxon, sci_name, specimen)
  
}

#' Tidy taxonomic data of pteridophytes of Japan
#'
#' Data is from Japan Green list
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
tidy_japan_names <- function (data) {
  data %>%
  select(taxon_id = ID20160331, scientific_name = `GreenList学名`,
         endemic = `固有`, conservation_status = `RL2012`) %>%
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
    left_join(japan_taxa) %>%
    # Convert full name with authors to only taxon, no authors
    assert(not_na, scientific_name) %>%
    left_join(
      parse_names_batch(.$scientific_name) %>% 
        select(scientific_name = b, taxon = c)
    ) %>%
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
#' This renames them to human-readable taxon names.
#' 
#' Uses parse_names_batch, which requires GNparser to be installed
#' and on PATH
#'
#' @param japan_rbcL rbcL alignment with names as codes
#' @param japan_taxa dataframe matching codes to species name
#'
rename_japan_rbcL_sexdip <- function (japan_rbcL, japan_taxa, repro_data) {
  
  # Make a table mapping dna taxon codes to scientific names
  # only 
  japan_sexdip_names_table <-
    tibble(
      original_name = names(japan_rbcL),
      taxon_id = names(japan_rbcL) %>%
        str_split("_") %>%
        map_chr(2)
    ) %>%
    left_join(japan_taxa) %>%
    # Convert full name with authors to only taxon, no authors
    assert(not_na, scientific_name) %>%
    left_join(
      parse_names_batch(.$scientific_name) %>% 
        select(scientific_name = b, taxon = c)
    ) %>%
    assert(not_na, taxon) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
    left_join(select(repro_data, taxon_id, sexual_diploid)) %>%
    filter(sexual_diploid == 1) %>%
    # Add JA so we know where it came from
    mutate(taxon = paste0(taxon, "_JAsexdip"))
  
  # Subset sequences to only sexual diploids
  japan_rbcL <- japan_rbcL[japan_sexdip_names_table$original_name]
  
  # Rename DNA sequences with taxon names
  names(japan_rbcL) <- japan_sexdip_names_table$taxon
  
  japan_rbcL
}

clean_taxonomy_data <- function (data) {
  
  data %>%
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
  
}

#' Read in a tree file contained in a zipped archive
#'
#' @param zip_folder Path to zip file
#' @param tree_file Name of nexus file within zip file
#'
#' @return List
#' 
read_tree_in_zip <- function (zip_folder, tree_file) {
  
  temp_dir <- tempdir()
  
  unzip(zip_folder, exdir = temp_dir)
  
  ape::read.tree(fs::path(temp_dir, tree_file))
  
}

#' Add new taxa to the Pteridophyte Phylogeny Group I (PPGI) data
#' 
#' PPGI was published in 2016. This adds new genera that have been
#' published since then.
#'
#' @param ppgi  Taxonomy of pteridophytes following Pteridophyte Phylogeny Group I
#'
#' @return Tibble
#'
add_new_pterido_taxa <- function (ppgi) {
  
  ppgi %>%
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
    )
  
}

# Checklist ----

make_checklist <- function (specimens, sci_names, taxonomy) {
  
  specimens %>%
    left_join(sci_names, by = "taxon") %>%
    left_join(taxonomy, by = "genus") %>%
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
    )
  
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
    filter(!is.na(date_collected)) %>%
    select(taxon, date_collected) %>%
    unique %>%
    mutate(abun = 1) %>%
    pivot_wider(names_from = date_collected, values_from = abun) %>%
    mutate_if(is.numeric, ~replace_na(., 0)) %>%
    column_to_rownames("taxon") %>%
    as.matrix %>%
    iNEXT(datatype = "incidence_raw", endpoint = endpoint)

  }

# Barcode analysis ----

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
#' @param rbcL_tree Phylogenetic tree
#' @param ppgi Taxonomy of pteridophytes following Pteridophyte Phylogeny Group I
#' @param outfile Name of file to use to save tree pdf
#'
#' @return Nothing; externally, the tree will be written as a pdf
#'
plot_rbcL_tree <- function(rbcL_tree, ppgi, specimens, dna_acc, outfile) {
  
  # Extract tips into tibble and add taxonomy
  tips <-
    tibble(tip = rbcL_tree$tip.label) %>%
    mutate(
      species = sp_name_only(tip, sep = "_"),
      genus = genus_name_only(tip, sep = "_")) %>%
    left_join(dplyr::select(ppgi, genus, family, class)) %>%
    mutate(genomicID = str_match(tip, "JNG.*$") %>% map_chr(1))
  
  # Identify lycophyte tips for rooting
  lycos <- tips %>% filter(class == "Lycopodiopsida") %>% pull(tip)
  
  # Root on lycophytes
  rbcL_tree <-
    ape::root(rbcL_tree, outgroup = lycos) %>%
    ape::ladderize(right = FALSE)
  
  # Reformat tip labels now that tip order has changed
  new_tips <-
    tibble(tip = rbcL_tree$tip.label) %>%
    mutate(genomicID = str_match(tip, "JNG.*$") %>% map_chr(1)) %>%
    left_join(dplyr::select(dna_acc, genomicID, specimen_id = specimenID)) %>%
    left_join(dplyr::select(specimens, specimen_id, taxon, specimen)) %>%
    mutate(new_tip = paste(taxon, specimen) %>% str_remove_all("Nitta ") %>% str_replace_all("_", " "))
  
  rbcL_tree$tip.label <- new_tips$new_tip
  
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