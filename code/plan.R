plan <- drake_plan(
  
  # Data loading and cleaning ----
  
  # Load PPGI taxonomy
  ppgi_path = target("data/ppgi_taxonomy.csv", format = "file"), # track file contents (not just path)
  ppgi = readr::read_csv(ppgi_path),
  
  # Load Nectandra specimen data
  nectandra_specimens_path = target("data/nectandra_specimens.csv", format = "file"),
  nectandra_specimens = readr::read_csv(nectandra_specimens_path),
  
  # Load Nectandra DNA accession data
  nectandra_dna_path = target("data/nectandra_DNA_accessions.csv", format = "file"),
  nectandra_dna = readr::read_csv(nectandra_dna_path),
  
  # Load Nectandra unaligned rbcL sequences
  nectandra_rbcL_raw_path = target("data/nectandra_rbcL.fasta", format = "file"),
  nectandra_rbcL_raw = ape::read.FASTA(nectandra_rbcL_raw_path),
  
  # Also read in one sequence that will be submitted to GenBank separately
  JNG4254_rbcL_raw_path = target("data/JNG4254.fasta", format = "file"),
  JNG4254_rbcL_raw = ape::read.FASTA(JNG4254_rbcL_raw_path),
  
  # Combine the unaligned Nectandra rbcL seqs
  nectandra_rbcL_with_JNG4254 = c(nectandra_rbcL_raw, JNG4254_rbcL_raw),
  
  # Load species richness and GPS locations of various protected
  # sites in Costa Rica
  cr_richness_path = target("data/costa_rica_richness.csv", format = "file"),
  cr_richness = readr::read_csv(cr_richness_path),
  
  # Unzip French Polynesia rbcL sequences
  # This requires doi_10.5061_dryad.df59g__v1.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g first
  moorea_rbcL_path = target({
    unzip_nitta_2017(
    dryad_zip_file = "data/doi_10.5061_dryad.df59g__v1.zip", 
    exdir = "data/nitta_2017")
    "data/nitta_2017/rbcL_clean_sporos.fasta"},
    format = "file"
  ),
  
  # Load French Polynesia rbcL sequences
  # (also add "_FP" to end of name)
  moorea_rbcL = ape::read.FASTA(moorea_rbcL_path) %>%
    purrr::set_names(., paste0(names(.), "_FP")),
  
  # Unzip Japan rbcL sequences and breeding mode data
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  ebihara_2019_data = target({
    unzip_ebihara_2019(
      dryad_zip_file = "data/doi_10.5061_dryad.4362p32__v4.zip", 
      exdir = "data/ebihara_2019"
    )
    c(tree = "data/ebihara_2019/rbcl_mrbayes.nex", 
      taxa = "data/ebihara_2019/FernGreenListV1.01E.xls",
      repro_data = "data/ebihara_2019/ESM1.csv")
  },
  format = "file"),
  
  # Load Japan rbcL data with names formatted as codes
  japan_rbcL_raw = ape::read.nexus.data(ebihara_2019_data[["tree"]]) %>% ape::as.DNAbin(),
  
  # Load table matching taxon codes to scientific names of Japanese pteridophytes
  japan_taxa = readxl::read_excel(ebihara_2019_data[["taxa"]]) %>% tidy_japan_names(),
  
  # Load data on reproductive mode for Japanese pteridophytes
  japan_repro_data = readr::read_csv(ebihara_2019_data[["repro_data"]]) %>%
    janitor::clean_names() %>%
    mutate(taxon_id = as.character(taxon_id)),
  
  # Rename Japan rbcL alignment as taxon names
  # (also add "_JA" to end of name)
  japan_rbcL = rename_japan_rbcL(
    japan_rbcL = japan_rbcL_raw, 
    japan_taxa = japan_taxa),
  
  # Also make an alignment of Japan sexual-diploids only
  japan_rbcL_sexdip = rename_japan_rbcL_sexdip(
    japan_rbcL = japan_rbcL_raw, 
    japan_taxa = japan_taxa, 
    japan_repro_data = japan_repro_data),
  
  # Checklist ----
  
  # Make species checklist, write out as SI
  checklist = make_checklist(
    specimens = nectandra_specimens, 
    taxonomy = ppgi) %>% 
    write_csv(file_out("results/table_S1.csv")),
  
  # Collection curve ----
  
  # Run iNEXT to generate interpolated/extrapolated species richness
  # using number of sampling days as the sampling unit
  # set endpoint (maximum number of collection days) to 150
  richness_estimate = estimate_richness_by_date(
    specimens = nectandra_specimens,
    endpoint = 150),
  
  # Barcode analysis ----
  
  # Align Nectandra sequences and add missing taxa
  nectandra_rbcL = align_rbcL(
    nectandra_rbcL_raw = nectandra_rbcL_with_JNG4254, 
    nectandra_dna = nectandra_dna, 
    nectandra_specimens = nectandra_specimens),
  
  # Calculate minimum interspecific distances for the three
  # rbcL datasets and bin them by 0.05% sequence divergence
  min_distance_table = analyze_min_dist(
    nectandra_rbcL = nectandra_rbcL, 
    moorea_rbcL = moorea_rbcL, 
    japan_rbcL = japan_rbcL, 
    japan_rbcL_sexdip = japan_rbcL_sexdip),
  
  # Phylogenetic analysis ----
  
  # Write out alignment for submission to Dryad
  rbcL_aln_out = phangorn::write.phyDat(
    x = nectandra_rbcL,
    file = file_out("results/nectandra_rbcL.phy"),
    format = "phylip"
  ),
  
  # Write out alignment for analysis with IQ-TREE
  rbcL_aln_out_iqtree = phangorn::write.phyDat(
    x = nectandra_rbcL, 
    file = file_out("iqtree_analysis/nectandra_rbcL"),
    format = "phylip"
  ),
  
  # Conduct ML phylogenetic analysis with IQ-TREE
  iqtree_results = jntools::iqtree(
    aln_path = file_in("iqtree_analysis/nectandra_rbcL"),
    wd = "iqtree_analysis",
    nt = 1,
    m = "TEST",
    bb = 1000,
    alrt = 1000,
    seed = 9130,
    redo = TRUE,
    echo = TRUE,
    tree_path = "iqtree_analysis/nectandra_rbcL.contree",
    produces1 = file_out("iqtree_analysis/nectandra_rbcL.treefile"),
    produced2 = file_out("iqtree_analysis/nectandra_rbcL.log")
  ),
  
  # Read in tree with SH-aLRT support (%) / UFboot support (%) at nodes.
  rbcL_tree = ape::read.tree(file_in("iqtree_analysis/nectandra_rbcL.treefile")),
  
  # Read in IQTREE log file to get stats about alignment and tree
  iqtree_log = readr::read_lines(file_in("iqtree_analysis/nectandra_rbcL.log")),
  
  # Print out tree for SI
  rbcL_tree_out = plot_rbcL_tree(
    phy = rbcL_tree,
    ppgi = ppgi,
    specimens = nectandra_specimens,
    dna_acc = nectandra_dna,
    outfile = file_out("results/Fig_S1.pdf")
  ),
  
  # Write out GenBank accession numbers for SI
  genbank_accession_table = make_genbank_accession_table(
    nectandra_rbcL = nectandra_rbcL, 
    DNA_accessions = nectandra_dna, 
    specimens = nectandra_specimens) %>% 
    write_csv(file_out("results/table_S2.csv")),
  
  # Render manuscript ----
  
  # First render to PDF, keeping the latex
  ms_pdf = render_tracked(
    knitr_in("ms/nectandra_pteridos.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results"),
    tracked_output = file_out(here::here("results/nectandra_pteridos.tex"))
  ),
  
  # Next use the latex to convert to docx with pandoc
  ms_docx = latex2docx(
    latex = file_in(here::here("results/nectandra_pteridos.tex")),
    docx = file_out(here::here("results/nectandra_pteridos.docx")),
    template = file_in(here::here("ms/plos-one.docx")),
    wd = here::here("results")
  ),
  
  # Also render the data readme for Dryad
  dryad_readme = render_tracked(
    knitr_in("ms/dryad_readme.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results"),
    tracked_output = file_out(here::here("results/dryad_readme.rtf"))
  ),
  
)
