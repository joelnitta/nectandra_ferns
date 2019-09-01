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

make_checklist <- function (specimens) {
  
  specimens %>%
    group_by(class, family, sci_name) %>%
    summarize(
      voucher = paste(specimen, collapse = ", ")
    )
  
}