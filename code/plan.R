plan <- drake_plan(
  
  # Unzip data ----
  
  # Unzip Nectandra ferns data
  # This requires doi_10.5061_dryad.bnzs7h477__v4.zip to be downloaded to data/
  # from https://doi.org/10.5061/dryad.bnzs7h477 first
  nec_ferns_data = target({
    unzip(
      zipfile = "data/doi_10.5061_dryad.bnzs7h477__v4.zip", 
      files = c(
        "costa_rica_richness.csv",
        "JNG4254.fasta",
        "nectandra_DNA_accessions.csv",
        "nectandra_gb_template.sbt",
        "nectandra_rbcL.fasta",
        "nectandra_specimens.csv",
        "ppgi_taxonomy.csv",
        "seqids.txt"),
      exdir = "data",
      overwrite = TRUE)
    c(
      costa_rica_richness.csv = "data/costa_rica_richness.csv",
      JNG4254.fasta = "data/JNG4254.fasta",
      nectandra_DNA_accessions.csv = "data/nectandra_DNA_accessions.csv",
      nectandra_gb_template.sbt = "data/nectandra_gb_template.sbt",
      nectandra_rbcL.fasta = "data/nectandra_rbcL.fasta",
      nectandra_specimens.csv = "data/nectandra_specimens.csv",
      ppgi_taxonomy.csv = "data/ppgi_taxonomy.csv",
      seqids.txt = "data/seqids.txt")
  },
  format = "file"),

  # Unzip Japan rbcL sequences and breeding mode data
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data/
  # from https://doi.org/10.5061/dryad.4362p32 first
  japan_ferns_data = target({
    unzip(
      zipfile = "data/doi_10.5061_dryad.4362p32__v4.zip", 
      files = c(
        "rbcl_mrbayes.nex",
        "FernGreenListV1.01E.xls",
        "ESM1.csv"),
      exdir = "data",
      overwrite = TRUE)
    c(
      rbcl_mrbayes.nex = "data/rbcl_mrbayes.nex",
      FernGreenListV1.01E.xls = "data/FernGreenListV1.01E.xls",
      ESM1.csv = "data/ESM1.csv")
  },
  format = "file"),
  
  # Unzip French Polynesia rbcL sequences
  # This requires doi_10.5061_dryad.df59g__v1.zip to be downloaded to data/
  # from https://doi.org/10.5061/dryad.df59g first
  moorea_rbcL_path = target({
    unzip_nitta_2017(
      dryad_zip_file = "data/doi_10.5061_dryad.df59g__v1.zip", 
      exdir = "data")
    "data/rbcL_clean_sporos.fasta"},
    format = "file"
  ),
  
  # Data loading and cleaning ----
  
  # Load PPGI taxonomy
  ppgi = readr::read_csv(nec_ferns_data[["ppgi_taxonomy.csv"]]),
  
  # Load Nectandra specimen data
  nectandra_specimens = readr::read_csv(nec_ferns_data[["nectandra_specimens.csv"]]),
  
  # Load Nectandra DNA accession data
  nectandra_dna = readr::read_csv(nec_ferns_data[["nectandra_DNA_accessions.csv"]]),
  
  # Load Nectandra unaligned rbcL sequences
  nectandra_rbcL_raw = ape::read.FASTA(nec_ferns_data[["nectandra_rbcL.fasta"]]),
  
  # Also read in one sequence that will be submitted to GenBank separately
  JNG4254_rbcL_raw = ape::read.FASTA(nec_ferns_data[["JNG4254.fasta"]]),
  
  # Load species richness and GPS locations of various protected
  # sites in Costa Rica
  cr_richness = readr::read_csv(nec_ferns_data[["costa_rica_richness.csv"]]),
  
  # Load Japan rbcL data with names formatted as codes
  japan_rbcL_raw = ape::read.nexus.data(japan_ferns_data[["rbcl_mrbayes.nex"]]) %>% ape::as.DNAbin(),
  
  # Load table matching taxon codes to scientific names of Japanese pteridophytes
  japan_taxa = readxl::read_excel(japan_ferns_data[["FernGreenListV1.01E.xls"]]) %>% tidy_japan_names(),
  
  # Load data on reproductive mode for Japanese pteridophytes
  japan_repro_data = readr::read_csv(japan_ferns_data[["ESM1.csv"]]) %>%
    janitor::clean_names() %>%
    mutate(taxon_id = as.character(taxon_id)),
  
  # Rename Japan rbcL alignment as taxon names
  # (also add "_JA" to end of name)
  japan_rbcL = rename_japan_rbcL(
    japan_rbcL = japan_rbcL_raw, 
    japan_taxa = japan_taxa),
  
  # Also make an alignment of Japan sexual-diploids only
  # (also add "_JAsexdip" to end of name)
  japan_rbcL_sexdip = rename_japan_rbcL_sexdip(
    japan_rbcL = japan_rbcL_raw, 
    japan_taxa = japan_taxa, 
    japan_repro_data = japan_repro_data),
  
  # Load French Polynesia rbcL sequences
  # (also add "_FP" to end of name)
  moorea_rbcL = ape::read.FASTA(moorea_rbcL_path) %>%
    purrr::set_names(., paste0(names(.), "_FP")),
  
  # Download elevation data for Costa Rica (for map)
  costa_rica_el = download_country_el("CRI"),
  
  # Checklist ----
  
  # Make species checklist, write out as manuscript SI
  checklist = make_checklist(
    specimens = nectandra_specimens, 
    taxonomy = ppgi) %>% 
    write_csv(s_table_path("checklist", ".csv")),
  
  # Format checklist for fernsoftheworld.com, write out
  fow_list = checklist %>%
    transmute(Family = family, Species = scientific_name, Collection = voucher, `Det.` = "J. Nitta") %>%
    separate_rows(Collection, sep = ", "),
  
  fow_list_out = write_csv(fow_list, "results/other/fow_nectandra_list.csv"),
  
  # Format collection data for fernsoftheworld.com, write out
  # - load list of accepted FOW genera
  fow_genera_path = target("data/fow_genera.txt", format = "file"),
  
  fow_genera = read_lines(fow_genera_path) %>%
    tibble(genus_with_auth = .) %>%
    mutate(genus = str_split(genus_with_auth, " ") %>% map_chr(1)),
  
  # - load list of accepted FOW species
  fow_species_path = target("data/fow_species.txt", format = "file"),
  
  fow_species = read_lines(fow_species_path) %>%
    tibble(species_with_auth = .) %>%
    taxastand::add_parsed_names(species_with_auth, species),
  
  fow_coll_data = format_coll_for_fow(nectandra_specimens, ppgi, fow_genera, fow_species),
  
  # Write out data for bulk upload to FOW blog posts one year at time
  fow_coll_data_out = fow_coll_data %>%
    mutate(year = str_split(fow_date, "-") %>% map_chr(3)) %>%
    assert(not_na, year) %>%
    verify(unique(year) == c(2008, 2011, 2013)) %>%
    group_split(year, keep = FALSE) %>%
    walk2(
      glue::glue("results/other/fow_nectandra_collection_data_{c(2008, 2011, 2013)}.csv"),
      ~write_csv(.x, .y, na = "")),
  
  # Collection curve ----
  
  # Run iNEXT to generate interpolated/extrapolated species richness
  # using number of sampling days as the sampling unit
  # set endpoint (maximum number of collection days) to 150
  richness_estimate = estimate_richness_by_date(
    specimens = nectandra_specimens,
    endpoint = 150),
  
  # Nectandra rbcL alignment ----
  
  # Get list of missing Nectandra taxa
  # (taxa that couldn't be successfully sequenced)
  nectandra_rbcL_missing_taxa = make_missing_taxa_list(
    nectandra_rbcL_raw = nectandra_rbcL_raw, 
    nectandra_dna = nectandra_dna, 
    nectandra_specimens = nectandra_specimens),
  
  # Download GenBank sequences for selected missing taxa
  genbank_rbcL_raw = fetch_gb_seqs(
    nectandra_rbcL_missing_taxa = nectandra_rbcL_missing_taxa,
    acc_keep = c(
      "KM008147", # Pteris altissima
      "AY175795", # Trichomanes polypodioides
      "U21289",   # Radiovittaria remota
      "AY095108"  # Abrodictyum ridigum
      )  
  ),
  
  # For now combine JNG4254 with other Nectandra seqs, but
  # this will need to be shifted to GenBank seqs when JNG4254 
  # accession number is ready (it was sequenced separately in the SIBN project)
  nectandra_rbcL_raw_with_JNG4254 = c(nectandra_rbcL_raw, JNG4254_rbcL_raw),
  
  # Rename newly sequenced Nectandra sequences
  nectandra_rbcL_raw_renamed = rename_nectandra_rbcL(
    nectandra_rbcL_raw = nectandra_rbcL_raw_with_JNG4254, 
    nectandra_dna = nectandra_dna,
    nectandra_specimens = nectandra_specimens
  ),
  
  # Combine and align newly sequenced and GenBank Nectandra sequences
  nectandra_rbcL = align_nectandra_rbcL(
    nectandra_rbcL_raw = nectandra_rbcL_raw_renamed, 
    genbank_rbcL_raw = genbank_rbcL_raw$seqs),
  
  # Barcode analysis ----
  
  # Remove Macrothelypteris torresiana (non-native species) from Nectandra alignment
  nectandra_rbcL_native_only = nectandra_rbcL[nectandra_rbcL %>% rownames %>% str_detect("Macrothelypteris_torresiana", negate = TRUE), ],
  
  # Combine the four rbcL datasets into a single alignment.
  # Differentiate datasets by code appended to each species name.
  combined_rbcL = align_all_rbcL(
    nectandra_rbcL = nectandra_rbcL_native_only, # _CR
    moorea_rbcL = moorea_rbcL, # _FP
    japan_rbcL = japan_rbcL, # _JA
    japan_rbcL_sexdip = japan_rbcL_sexdip # _JAsexdip
  ),
  
  # Calculate minimum interspecific distances for each
  # dataset separately, bin them by 0.05% sequence divergence,
  # and combine these into a single dataframe
  min_distance_table = purrr::map_df(
    c("CR", "FP", "JA", "JAsexdip"), 
    ~bin_min_inter_dist_by_dataset(
      full_aln = combined_rbcL, 
      dataset_select = .)
  ),
  
  # Nectandra rbcL tree ----
  
  # Write out rbcL alignment for submission to Dryad
  rbcL_aln_out = phangorn::write.phyDat(
    x = nectandra_rbcL,
    file = file_out("results/dryad/nectandra_rbcL.phy"),
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
  
  # Print out tree for manuscript SI
  rbcL_tree_out = plot_nectandra_rbcL_tree(
    phy = rbcL_tree,
    ppgi = ppgi,
    outfile = s_fig_path("tree", ".pdf")
  ),
  
  # Load new GenBank accession numbers
  # (MW138110 - MW138295, plus MT657442 submitted separately)
  new_genbank_accs = load_genbank_accs(nec_ferns_data[["seqids.txt"]]),
  
  # Write out GenBank accession numbers table for manuscript SI
  genbank_accession_table = make_genbank_accession_table(
    new_genbank_accs = new_genbank_accs,
    nectandra_rbcL_raw = nectandra_rbcL_raw_with_JNG4254, 
    DNA_accessions = nectandra_dna, 
    specimens = nectandra_specimens,
    genbank_rbcL_metadata = genbank_rbcL_raw$metadata,
    nectandra_rbcL_aligned = nectandra_rbcL) %>% 
    write_csv(s_table_path("genbank", ".csv")),
  
  # Other rbcL trees----
  
  # - Cyatheaceae: download sequences from GenBank and make alignment, write out for Dryad
  cyatheaceae_seqs = make_broad_alignment(
    nectandra_rbcL = nectandra_rbcL, 
    ppgi = ppgi,
    nectandra_family = c("Cyatheaceae", "Dicksoniaceae"),  # Include Dicksoniaceae as outgroup
    genbank_group = "Cyatheaceae"
  ),
  
  cyatheaceae_seqs_out = ape::write.FASTA(cyatheaceae_seqs, file_out("results/dryad/cyatheaceae_rbcL.fasta")),
  
  # - Cyatheaceae: Make tree with fasttree, write out for Dryad
  cyatheaceae_tree = fasttree(cyatheaceae_seqs),
  
  cyatheaceae_tree_out = ape::write.tree(cyatheaceae_tree, file_out("results/dryad/cyatheaceae_rbcL.tre")),
  
  # - Cyatheaceae: Identify outgroup sequences
  cyatheaceae_outgroup = identify_outgroup(
    phy = cyatheaceae_tree, 
    ppgi = ppgi,
    family == "Dicksoniaceae"
  ),
  
  # - Cyatheaceae: Print out tree for manuscript SI
  cyatheaceae_rbcL_tree_pdf = plot_broad_rbcL_tree(
    phy = cyatheaceae_tree,
    outgroup = cyatheaceae_outgroup,
    outfile = s_fig_path("cyatheaceae-tree", ".pdf")
  ),
  
  # - Grammitids: download sequences from GenBank and make alignment, write out for Dryad
  grammitid_seqs = make_broad_alignment(
    nectandra_rbcL = nectandra_rbcL, 
    ppgi = ppgi,
    nectandra_family = "Polypodiaceae", # Include other Polypodiaceae as outgroup
    genbank_group = "Grammitidoideae",
    exclude_list = "MH159215" # Exclude misidentified Ascogrammitis anfractuosa on GenBank
  ),
  
  grammitid_seqs_out = ape::write.FASTA(grammitid_seqs, file_out("results/dryad/grammitidoideae_rbcL.fasta")),
  
  # - Grammitids: Make tree with fasttree, write out for Dryad
  grammitid_tree = fasttree(grammitid_seqs),
  
  grammitid_tree_out = ape::write.tree(cyatheaceae_tree, file_out("results/dryad/grammitidoideae_rbcL.tre")),
  
  # - Grammitids: Identify outgroup sequences
  grammitid_outgroup = identify_outgroup(
    phy = grammitid_tree, 
    ppgi = ppgi,
    subfamily == "Polypodioideae"
  ),
  
  # - Grammitids: Print out tree for manuscript SI
  grammitid_rbcL_tree_pdf = plot_broad_rbcL_tree(
    phy = grammitid_tree,
    outgroup = grammitid_outgroup,
    nodelab_size = 0.8,
    outfile = s_fig_path("grammitid-tree", ".pdf")
  ),
  
  # Manuscript rendering ----
  
  # Track bibliography files
  refs = target("ms/references.bib", format = "file"),
  refs_other = target("ms/references_other.yaml", format = "file"),
  
  # First render to PDF, keeping the latex
  ms_pdf = render_tracked(
    input = knitr_in("ms/nectandra_ferns.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results/ms"),
    tracked_output = file_out(here::here("results/ms/nectandra_ferns.tex")),
    dep1 = refs,
    dep2 = refs_other
  ),
  
  # Next use the latex to convert to docx with pandoc
  ms_docx = latex2docx(
    latex = file_in(here::here("results/ms/nectandra_ferns.tex")),
    docx = file_out(here::here("results/ms/nectandra_ferns.docx")),
    template = file_in(here::here("ms/plos-one.docx")),
    wd = here::here("results/ms")
  ),
  
  # Also render the data readme for Dryad
  dryad_readme = render_tracked(
    input = knitr_in("ms/dryad_readme.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results/dryad"),
    tracked_output = file_out(here::here("results/dryad/dryad_readme.rtf"))
  ),
  
  # GenBank submission ----
  
  # Format genbank metadata
  genbank_metadata = format_data_for_genbank(
    sample_data = nectandra_specimens, 
    dna_data = nectandra_dna, 
    seqs = nectandra_rbcL_raw, 
    adjust_remainder = c("0" = 2, "1" = 1, "2" = 3), 
    manual_readframe = NULL, 
    notes = NULL
  ),
  
  # Format fasta headers for tbl2asn.
  # These contain all the necessary metadata for each sample
  seqs_for_tbl2asn = format_fasta_for_tbl2asn(genbank_metadata, nectandra_rbcL_raw),
  
  # Make feature table (describing genes) to use for tbl2asn
  genbank_features = genbank_metadata %>%
    mutate(feature = make_feature(name = genomic_id, seq_len = seq_len, codon_start = codon_start)) %>%
    pull(feature),
  
  # Run tbl2asn
  tbl2asn_results = tbl2asn(
    genbank_template_file = nec_ferns_data[["nectandra_gb_template.sbt"]],
    genbank_features = genbank_features, 
    seqs = seqs_for_tbl2asn, 
    submission_name = "nectandra_ferns_rbcL", 
    results_dir = "results/genbank_submission"),
  
  # Check that amino acid translations are correct
  
  # (AA seqs in GenBank submission are based on the assumption that all DNA sequences 
  # are exactly the correct length, and translations made by shifting reading frame appropriately)
  flatfile_path = target(
    tracked_file("results/genbank_submission/nectandra_ferns_rbcL.gbf", tbl2asn_results),
    format = "file"),
  
  # - parse translated amino acids from flatfile
  aa_seqs = parse_aa_from_flatfile(flatfile_path)
  
)
