plan <- drake_plan(
  
  ppgi = read_csv("data_raw/ppgi_taxonomy.csv"),
  
  taxonomy = tidy_taxonomy(ferncollectR::taxonomy_data$species),

  specimens = tidy_specimens(ferncollectR::collection_data$specimens, ppgi, taxonomy),
    
  checklist = make_checklist(specimens)
    
)