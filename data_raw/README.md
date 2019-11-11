# Data_raw

Raw data files used in analysis.

costa_rica_richness.csv: Data on species richness of pteridophytes in protected 
areas in Costa Rica. Compiled by hand by Joel Nitta based on references in the 
"citation" column. Numbers for La Selva calculated using script in 
code/count_la_selva_taxa.R.

DNA_accessions.csv: DNA accessions by Joel Nitta, including Nectandra and other 
sites. This is filtered to create data/nectrandra_dna.csv. Copied from 
/Users/joelnitta/Google Drive/collections_db.

doi_10.5061_dryad.df59g__v1.zip: Zip file of data from Nitta, Joel H.; Meyer, 
Jean-Yves; Taputuarai, Ravahere; Davis, Charles C. (2016) Life cycle matters: 
DNA barcoding reveals contrasting community structure between fern sporophytes 
and gametophytes, Dryad, Dataset, https://doi.org/10.5061/dryad.df59g

doi_10.5061_dryad.4362p32__v4.zip: Zip file of data from Ebihara, Atsushi; 
Nitta, Joel H. (2019) An update and reassessment of fern and lycophyte diversity 
data in the Japanese Archipelago, v5, Dryad, Dataset, 
https://doi.org/10.5061/dryad.4362p32

geneious_sporos_rbcL.fasta: Consensus rbcL sequences for ferns and lycophytes
of Nectandra (sporophytes). Sequences were exported from Geneious as follows:
  1. select all clean sequences in a folder (here, Nectandra Clean Sporos)
  2. File -> Export -> selected documents -> as fasta -> 
  save as geneious_sporos_rbcL.fasta in data-raw
Note that each sequence will be named using the "name" column in Genenious. 
This should just be my unique DNA sequence ID number ("g_number")

Lista_especies_LS_feb2017.xlsx: Checklist of all plant species at the La Selva 
protected area in Costa Rica. Downloaded from https://sura.ots.ac.cr/florula4/ 
on 2019-09-03.

photos.csv: Plant voucher photo file names by Joel Nitta, including Nectandra 
and other sites. Copied from /Users/joelnitta/Google Drive/collections_db.

ppgi_taxonomy.csv: Spreadsheet of the Pteridophyte Phylogeny I working group
taxonomic system for pteridophytes at the genus level and above 
(The Pteridophyte Phylogeny Group, 2016. A community-derived classification for 
extant lycophytes and ferns. J Syst Evol 54:563-606). Includes columns for 
class, order, suborder, family, subfamily, and genus. Encoding is Unicode (UTF-8).

specimens.csv: Specimen collection data of all collections by Joel Nitta, 
including Nectandra and other sites. This is filtered to create 
data/nectrandra_specimens.csv. 
Copied from /Users/joelnitta/Google Drive/collections_db.

taxonomy.csv: Taxonomic data for all collections by Joel Nitta. from Nectandra 
and other sites.  Copied from /Users/joelnitta/Google Drive/collections_db.
