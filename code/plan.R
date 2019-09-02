plan <- drake_plan(
  
  ppgi = read_csv("data_raw/ppgi_taxonomy.csv"),
  
  taxonomy = tidy_taxonomy(ferncollectR::taxonomy_data$species),

  specimens = tidy_specimens(ferncollectR::collection_data$specimens, ppgi, taxonomy),
    
  checklist = make_checklist(specimens),
  
  # Run iNEXT to generate interpolated/extrapolated species richness
  # set endpoint (maximum number of individuals projected to collected) to 800
  richness_estimate = count(specimens, taxon) %>%
    pull(n) %>%
    iNEXT(datatype="abundance", endpoint=800)
  
)