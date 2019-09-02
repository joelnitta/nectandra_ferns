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
  ggplot (inext_out$iNextEst, aes(x=m, y=qD)) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), fill = "grey70", alpha=0.2) +
    geom_line(aes(linetype=method), size=0.8) +
    scale_linetype_manual(values=c("dotted", "solid", "solid")) +
    geom_point(data=observed_data, size=2.5) +
    scale_y_continuous("Richness") +
    scale_x_continuous("Number of individuals") +
    standard_theme2() +
    theme(legend.position="none", 
          legend.title=element_blank())
}
