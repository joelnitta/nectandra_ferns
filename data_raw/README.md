# Data

Data files used in analysis.

DNA_accessions.csv: DNA accessions by Joel Nitta, including Nectandra and other 
 sites. This is filtered to create data/nectrandra_dna.csv. Copied from 
/Users/joelnitta/Google Drive/collections_db.

ESM1.csv: A list of native fern and lycophyte taxa (species, subspecies and varieties) in Japan. Taxon ID refers to that in FernGreenList ver.1.0.1 (http://www.rdplants.org/gl/). Unless otherwise noted, rbcL GenBank accession numbers 
are those used in Ebihara et al. (2010). Asterisks indicate newly generated sequences 
by this study. Voucher information only provided for newly generated sequences. 
Information on reproductive modes and ploidy levels follow those in Ebihara et 
al. (2016, 2017), and only records based on material collected in Japan are used. 
For reproductive mode, irregular meiosis is not considered, 0 = no information, 
1 = sexual, 2 = apomictic and 3 = sexual + apomictic. Encoding is Unicode (UTF-8).

FernGreenListV1.01.xls: List of Japanese ferns and lycophytes species including
scientific name, endemic status, conservation status, and other taxonomic data.
Downloaded from http://www.rdplants.org/gl/FernGreenListV1.01.xls on 2019-07-17.

ppgi_taxonomy.csv: Spreadsheet of the Pteridophyte Phylogeny I working group
taxonomic system for pteridophytes at the genus level and above (The Pteridophyte
Phylogeny Group, 2016. A community-derived classification for extant lycophytes
and ferns. J Syst Evol 54:563-606). Includes columns for class, order, suborder,
family, subfamily, and genus. Encoding is Unicode (UTF-8).

rbcl_mrbayes.nex: NEXUS file used for phylogenetic analysis of Japanese fern
and lycophyte taxa with MrBayes. From Ebihara and Nitta 2019 J Plant Res.

specimens.csv: Specimen collection data of all collections by Joel Nitta, 
including Nectandra and other sites. This is filtered to create 
data/nectrandra_specimens.csv. 
Copied from /Users/joelnitta/Google Drive/collections_db.

sporos_rbcL.phy: Cleaned sporophyte rbcL sequences. Copied from R/nectandraferns/data_raw.

taxonomy.csv: Taxonomic data for all collections by Joel Nitta, including taxa 
from Nectandra and other sites. 
Copied from /Users/joelnitta/Google Drive/collections_db.

geneious_sporos_rbcL.fasta: Consensus rbcL sequences for ferns and lycophytes
of Nectandra (sporophytes). Sequences were exported from Geneious as follows:
  1. select all clean sequences in a folder (here, Nectandra Clean Sporos)
  2. File -> Export -> selected documents -> as fasta -> save as geneious_sporos_rbcL.fasta in data-raw
Note that each sequence will be named using the "name" column in Genenious. This should just be my
unique DNA sequence ID number ("g_number")

geneious_seq_quality.csv: Sequence quality data for geneious_sporos_rbcL.fasta
exported from Geneious.
Sequence quality metadata was exported from Geneious as follows:
  1. select all clean sequences in a folder (here, Nectandra Clean Sporos)
  2. File -> Export -> selected documents -> as csv -> select the following columns:
  Name, Ambiguities, HQ%, MQ%, LQ%, Sequence length -> save as geneious_seq_quality.csv in data-raw
This could also include sequence quality from other projects as long as Name (G#) is unique.
Ambiguities is number of ambiguous bases (non ATCG, I think) in the sequence
% columns are percentage of bases with quality score in that bin (high, med, low quality).
