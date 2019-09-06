plan <- drake_plan(
  
  # Data loading and cleaning ----
  
  # Load PPGI taxonomy
  ppgi = read_csv("data_raw/ppgi_taxonomy.csv"),
  
  # Load Nectandara taxonomy data with scientific names including authors
  sci_names = read_csv(file_in("data_raw/taxonomy.csv")) %>% 
    clean_taxonomy_data,

  # Load pre-processed Nectandra specimen data
  specimens = read_csv(file_in("data/nectandra_specimens.csv")),
  
  # Load DNA accession data
  dna_acc = read_csv("data_raw/DNA_accessions.csv"),
  
  # Download French Polynesia rbcL sequences
  nitta_2017_data = download_and_unzip_nitta_2017(
    dl_path = "data_raw/nitta_2017_data_and_scripts.zip",
    unzip_path = "data_raw/nitta_2017/",
    out1 = file_out("data_raw/nitta_2017/rbcL_clean_sporos.fasta")),
  
  # Load French Polynesia rbcL sequences
  # (also add "_FP" to end of name)
  moorea_rbcL = read.FASTA(file_in("data_raw/nitta_2017/rbcL_clean_sporos.fasta")) %>%
    purrr::set_names(., paste0(names(.), "_FP")),
  
  # Load Japan rbcL data with names formatted as codes
  japan_rbcL_raw = read.nexus.data("data_raw/rbcl_mrbayes.nex") %>% as.DNAbin,
  
  # Load table matching taxon codes to scientific names of Japanese pteridophytes
  japan_taxa = read_excel("data_raw/FernGreenListV1.01.xls") %>% tidy_japan_names(),
  
  # Load data on reproductive mode for Japanese pteridophytes
  repro_data = read_csv("data_raw/ESM1.csv") %>%
    clean_names %>%
    mutate(taxon_id = as.character(taxon_id)),
  
  # Rename Japan rbcL alignment as taxon names
  # (also add "_JA" to end of name)
  japan_rbcL = rename_japan_rbcL(japan_rbcL_raw, japan_taxa),
  
  # Also make an alignment of Japan sexual-diploids only
  japan_rbcL_sexdip = rename_japan_rbcL_sexdip(
    japan_rbcL_raw, japan_taxa, repro_data),
  
  # Load Nectandra rbcL sequences
  # (also add "_CR" to end of name)
  nectandra_rbcL = read.dna("data/nectandra_rbcL.phy") %>% 
    as.list %>%
    set_names(., paste0(names(.), "_CR")),
  
  # Load species richness and GPS locations of various protected
  # sites in Costa Rica
  cr_richness = read_csv("data_raw/costa_rica_richness.csv"),
  
  # Checklist ----
  # Make species checklist, write out as SI
  checklist = make_checklist(specimens, sci_names, ppgi) %>% 
    write_csv(here("ms/table_S1.csv")),
  
  # Collection curve ----
  # Run iNEXT to generate interpolated/extrapolated species richness
  # set endpoint (maximum number of individuals projected to collected) to 1000
  richness_estimate = count(specimens, taxon) %>%
    pull(n) %>%
    iNEXT(datatype="abundance", endpoint=1000),
  
  # Barcode analysis ----
  
  # Calculate minimum interspecific distances for the three
  # rbcL datasets and bin them by 0.05% sequence divergence
  min_distance_table = analyze_min_dist(
    nectandra_rbcL, moorea_rbcL, japan_rbcL, japan_rbcL_sexdip),
  
  # Phylogenetic analysis ----
  
  # Load Nectandra rbcL ML tree
  # (output of running RAxML on rbcL alignment on CIPRES)
  rbcL_tree = read_tree_in_zip("data/nectandra_rbcL_cipres.zip", "RAxML_bipartitions.result"),
  
  # Print out tree for SI
  rbcL_tree_out = plot_rbcL_tree(
    rbcL_tree,
    ppgi,
    specimens,
    dna_acc,
    file_out(here("ms/Fig_S1.pdf"))
  ),
  
  # Render manuscript ----
  ms = rmarkdown::render(
    knitr_in(here("ms/nectandra_pteridos.Rmd")),
    output_file = file_out(here("ms/nectandra_pteridos.pdf")),
    quiet = TRUE)
  
)
