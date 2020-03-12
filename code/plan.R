plan <- drake_plan(
  
  # Data loading and cleaning ----
  
  # Load PPGI taxonomy
  ppgi = read_csv(file_in("data/ppgi_taxonomy.csv")),
  
  # Load Nectandra specimen data
  specimens = read_csv(file_in("data/nectandra_specimens.csv")),
  
  # Load DNA accession data
  DNA_accessions = read_csv(file_in("data/nectandra_DNA_accessions.csv")),
  
  # Load unaligned rbcL sequences
  nectandra_rbcL_raw = read.FASTA("data/nectandra_rbcL.fasta"),
  
  # Load species richness and GPS locations of various protected
  # sites in Costa Rica
  cr_richness = read_csv(file_in("data/costa_rica_richness.csv")),
  
  # Unzip French Polynesia rbcL sequences
  # This requires doi_10.5061_dryad.df59g__v1.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g first
  nitta_2017_data = unzip_nitta_2017(
    dryad_zip_file = file_in("data/doi_10.5061_dryad.df59g__v1.zip"), 
    exdir = "data/nitta_2017",
    produces = file_out("data/nitta_2017/rbcL_clean_sporos.fasta")
    ),
  
  # Load French Polynesia rbcL sequences
  # (also add "_FP" to end of name)
  moorea_rbcL = read.FASTA(file_in("data/nitta_2017/rbcL_clean_sporos.fasta")) %>%
    purrr::set_names(., paste0(names(.), "_FP")),
  
  # Unzip Japan rbcL sequences and breeding mode data
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  ebihara_2019_data = unzip_ebihara_2019(
    dryad_zip_file = file_in("data/doi_10.5061_dryad.4362p32__v4.zip"), 
    exdir = "data/ebihara_2019",
    produces_1 = file_out("data/ebihara_2019/rbcl_mrbayes.nex"),
    produces_2 = file_out("data/ebihara_2019/FernGreenListV1.01E.xls"),
    produces_3 = file_out("data/ebihara_2019/ESM1.csv")
  ),
  
  # Load Japan rbcL data with names formatted as codes
  japan_rbcL_raw = read.nexus.data(file_in("data/ebihara_2019/rbcl_mrbayes.nex")) %>% as.DNAbin,
  
  # Load table matching taxon codes to scientific names of Japanese pteridophytes
  japan_taxa = read_excel(file_in("data/ebihara_2019/FernGreenListV1.01E.xls")) %>% tidy_japan_names(),
  
  # Load data on reproductive mode for Japanese pteridophytes
  repro_data = read_csv(file_in("data/ebihara_2019/ESM1.csv")) %>%
    clean_names %>%
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
    repro_data = repro_data),

  # Checklist ----
  
  # Make species checklist, write out as SI
  checklist = make_checklist(specimens, ppgi) %>% 
    write_csv(file_out("ms/table_S1.csv")),
  
  # Collection curve ----
  
  # Run iNEXT to generate interpolated/extrapolated species richness
  # using number of sampling days as the sampling unit
  # set endpoint (maximum number of collection days) to 150
  richness_estimate = estimate_richness_by_date(
    specimens,
    endpoint = 150),
  
  # Barcode analysis ----
  
  # Align Nectandra sequences and add missing taxa
  nectandra_rbcL = align_rbcL(
    nectandra_rbcL_raw, 
    DNA_accessions, 
    specimens),
  
  # Calculate minimum interspecific distances for the three
  # rbcL datasets and bin them by 0.05% sequence divergence
  min_distance_table = analyze_min_dist(
    nectandra_rbcL, moorea_rbcL, japan_rbcL, japan_rbcL_sexdip),
  
  # Phylogenetic analysis ----
  
  # Conduct ML phylogenetic analysis with IQ-TREE
  iqtree_results = jntools::iqtree(
    alignment = nectandra_rbcL,
    wd = "iqtree_analysis",
    nt = 1,
    m = "TEST",
    bb = 1000,
    alrt = 1000,
    seed = 9130,
    redo = TRUE,
    echo = TRUE,
    produces1 = file_out("iqtree_analysis/nectandra_rbcL.phy.treefile"),
    produced2 = file_out("iqtree_analysis/nectandra_rbcL.phy.log")
  ),
  
  rbcL_tree = ape::read.tree(file_in("iqtree_analysis/nectandra_rbcL.phy.treefile")),
  
  iqtree_log = read_lines(file_in("iqtree_analysis/nectandra_rbcL.phy.log")),
  
  # Write out alignment for dryad
  rbcL_aln_out = phangorn::write.phyDat(nectandra_rbcL, "results/nectandra_rbcL.phy"),
  
  # Print out tree for SI
  rbcL_tree_out = plot_rbcL_tree(
    rbcL_tree,
    ppgi,
    specimens,
    DNA_accessions,
    file_out("ms/Fig_S1.pdf")
  ),
  
  # Write out GenBank accession numbers for SI
  genbank_accession_table = make_genbank_accession_table(
    nectandra_rbcL, DNA_accessions, specimens) %>% 
    write_csv(file_out("ms/table_S2.csv")),
  
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
  )
  
)
