# Set up captions

# - Figures
figure_full <- captioner::captioner(prefix = "Fig", suffix = ".")
figure_full(name = "map", caption = 'Location of the Nectandra Cloud Forest Reserve ("Nectandra").')
figure_full(name = "richness-cr", caption = "Species richness of pteridophytes at protected areas in Costa Rica.")
figure_full(name = "richness-inext", caption = "Interpolation (solid line) and extrapolation (dashed line) of species richness of pteridophytes at Nectandra Cloud Forest Reserve, Costa Rica.")
figure_full(name = "min-dist", caption = "Minimum interspecific distances for selected pteridophyte floras.")

# - Tables
# table_full <- captioner::captioner(prefix = "Table")

# - SI figures
s_figure_full <- captioner::captioner(prefix = "S", auto_space = FALSE, suffix = " Fig. ")
s_figure_full(name = "tree", "Maximum-likelihood phylogenetic tree of ferns and lycophytes at the Nectandra Cloud Forest Reserve, Costa Rica inferred using IQ-TREE with 1,000 SH-like approximate likelihood ratio test (SH-aLRT) and ultra-fast bootstrap (UFboot) replicates each.")
s_figure_full(name = "cyatheaceae-tree", "Maximum-likelihood phylogenetic tree of the family Cyatheaceae including all available *rbcL* sequences on GenBank and newly sequenced taxa from the Nectandra Cloud Forest Reserve, Costa Rica inferred using FastTree.")
s_figure_full(name = "grammitid-tree", "Maximum-likelihood phylogenetic tree of grammitid ferns (subfamily Grammitidoideae) including all available *rbcL* sequences on GenBank and newly sequenced taxa from the Nectandra Cloud Forest Reserve, Costa Rica inferred using FastTree.")

# - SI tables
s_table_full <- captioner::captioner(prefix = "S", auto_space = FALSE, suffix = " Table. ")
s_table_full(name = "genbank", "GenBank accession numbers of sequences analyzed in this study.")
s_table_full(name = "checklist", "A checklist of fern and lycophyte species observed at the Nectandra Cloud Forest Reserve, Costa Rica.")

# Make short versions of citation functions
# - Just the number
figure <- pryr::partial(figure_full, display = "cite")
# table <- pryr::partial(table_full, display = "cite")
s_figure <- function(x) {s_figure_full(x) %>% stringr::str_match("(S[0-9] Fig)") %>% magrittr::extract(,2)}
s_table <- function(x) {s_table_full(x) %>% stringr::str_match("(S[0-9] Table)") %>% magrittr::extract(,2)}

# - Just the caption
figure_cap <- function(x) {figure_full(x) %>% stringr::str_match("Fig [0-9]+\\. (.*)$") %>% magrittr::extract(,2)}
# table_cap <- function(x) {table_full(x) %>% stringr::str_match("Table [0-9]+ (.*)$") %>% magrittr::extract(,2)}

# Optional:
# Check the order of figure and table citations in the Rmd file.
# Adjust the order of creation of captions above if it doesn't match.

# check_citation_order("ms/nectandra_pteridos.Rmd", "figure")

# check_citation_order("ms/nectandra_pteridos.Rmd", "s_figure")
 
# check_citation_order("ms/nectandra_pteridos.Rmd", "s_table")
