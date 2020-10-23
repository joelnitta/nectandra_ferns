# Data_raw

Raw data files not uploaded to Dryad.

To copy DNA_accessions.csv, photos.csv, specimens.csv, and taxonomy.csv into this folder, run:
cp -p /Users/joelnitta/Google\ Drive/collections_db/*.csv /Users/joelnitta/repos/nectandra/data_raw/

costa_rica_richness_raw.csv: Data on species richness of pteridophytes in protected 
areas in Costa Rica. Compiled by hand by Joel Nitta based on references in the 
"citation" column. Does not include richness for La Selva, which is calculated from
Lista_especies_LS_feb2017.xlsx.

DNA_accessions.csv: DNA accessions by Joel Nitta, including Nectandra and other 
sites. This is filtered to create data/nectandra_DNA_accessions.csv. Copied from 
/Users/joelnitta/Google Drive/collections_db.

Lista_especies_LS_feb2017.xlsx: Checklist of all plant species at the La Selva 
protected area in Costa Rica. Downloaded from https://sura.ots.ac.cr/florula4/ 
on 2019-09-03.

photos.csv: Plant voucher photo file names by Joel Nitta, including Nectandra 
and other sites. Copied from /Users/joelnitta/Google Drive/collections_db.

ppgi_taxonomy_raw.csv: Spreadsheet of the Pteridophyte Phylogeny I working group
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
