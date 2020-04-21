# Data

Data files used in this analysis.

--------------------------

costa_rica_richness.csv: Data on species richness of pteridophytes in protected
areas in Costa Rica. Compiled by Joel Nitta based on references in the
“citation” column.

Number of variables: 12

Number of cases/rows: 6

Variable list:

•	name: Abbreviated name of site.

•	full_name: Full name of site.

•	min_el_m: Minimum elevation of site in meters.

•	max_el_m: Maximum elevation of site in meters.

•	area_ha: Area of site in hectares.

•	richness: Number of species occurring at the site.

•	richness_per_ha: Number of species per hectare occurring at the site.

•	holdridge_type: Holdridge (1967) life-zone type.

•	citation: Reference for data.

•	citation_number: Reference number in manuscript.

•	latitude: Latitude in decimal-degrees.

•	longitude: Longitude in decimal-degrees.

Missing data codes: Missing data have no values (nothing entered between commas
in the CSV file).

Specialized formats or other abbreviations used: none.

--------------------------

doi_10.5061_dryad.df59g__v1.zip: Zip file of data from Nitta, Joel H.; Meyer,
Jean-Yves; Taputuarai, Ravahere; Davis, Charles C. (2016) Life cycle matters:
DNA barcoding reveals contrasting community structure between fern sporophytes
and gametophytes, Dryad, Dataset, https://doi.org/10.5061/dryad.df59g

--------------------------

doi_10.5061_dryad.4362p32__v4.zip: Zip file of data from Ebihara, Atsushi;
Nitta, Joel H. (2019) An update and reassessment of fern and lycophyte diversity
data in the Japanese Archipelago, v5, Dryad, Dataset,
https://doi.org/10.5061/dryad.4362p32

--------------------------

JNG4254.fasta: DNA sequence in FASTA format of rbcL gene from Amauropelta
atrovirens (C. Chr.) Salino & T.E. Almeida (Nitta 2237).

--------------------------

nectandra_DNA_accessions.csv: DNA accession numbers and specimen accession
numbers of pteridophytes from the Nectandra Cloud Forest Reserve, Costa Rica.

Number of variables: 2

Number of cases/rows: 237

Variable list:

•	genomic_id: Genomic accession number assigned during DNA extraction, of the
form “JNG” plus a four-digit number. Unique values.

•	specimen_id: Specimen accession number assigned to each specimen in
nectandra_specimens.csv. Integer (not unique).

Missing data codes: no missing data.

Specialized formats or other abbreviations used: none.

--------------------------

nectandra_rbcL.fasta: Newly generated rbcL sequences of pteridophytes from the
Nectandra Cloud Forest Reserve, Costa Rica in FASTA format. All species included
occur at the Nectandra Cloud Forest Reserve, Costa Rica; a small number of
sequences are from specimens collected elsewhere. Sequence names correspond to
‘genomic_id’ in nectandra_DNA_accessions.csv. 188 sequences; shortest sequence
466 bp; longest sequence 1309 bp; mean sequence length 1291 bp. Exported from
Geneious project folder “Clean Sporos Trimmed Genbank Submission” (raw Geneious
project file not included in this dataset).

--------------------------

nectandra_specimens.csv: Specimen data of pteridophytes from the Nectandra Cloud
Forest Reserve, Costa Rica collected by Joel Nitta. Formatting UTF-8.

Number of variables: 21

Number of cases/rows: 322

Variable list:

•	specimen_id: Unique specimen identification number (integer).

•	specimen: Voucher specimen number.

•	genus: Genus

•	specific_epithet: Specific epithet.

•	infraspecific_rank: Infraspecific rank.

•	infraspecific_name: Infraspecific name.

•	certainty: Degree of taxonomic certainty if not completely certain.

•	species: Species (genus plus specific epithet).

•	taxon: Species plus infraspecific name.

•	scientific_name: Taxon plus scientific author.

•	country: Country of origin.

•	locality: General area of collection.

•	site: Specific site where collected.

•	observations: Observations about specimen.

•	elevation: Elevation in m.

•	latitude: Latitude in decimal-degrees.

•	longitude: Longitude in decimal-degrees.

•	collector: Name of collector.

•	other_collectors: Names of other collectors if present.

•	herbaria: Codes of herbaria where voucher specimens are lodged.

•	date_collected: Date collected in YYYY-MM-DD format.

Missing data codes: Missing or non-applicable data have no values (nothing
entered between commas in the CSV file).

Specialized formats or other abbreviations used: Herbaria codes follow Index
Herbariorum (Thiers 2020), except for “Nectandra”, which indicates the private
herbarium at the Nectandra Cloud Forest Reserve.

--------------------------

ppgi_taxonomy.csv: Taxonomic system of Pteridophyte Phylogeny Group I (2016) for
pteridophytes at the genus level and above. Updated with one new genus (Hiya).

Number of variables: 6

Number of cases/rows: 338

Variable list:

•	class: Class.

•	order: Order.

•	suborder: Suborder.

•	family: Family.

•	subfamily: Subfamily.

•	genus: Genus.

Missing data codes: Non-applicable data have no values (nothing entered between
commas in the CSV file).

Specialized formats or other abbreviations used: none.

--------------------------

REFERENCES

Pteridophyte Phylogeny Group I (2016) A community-derived classification for
extant lycophytes and ferns. Journal of Systematics and Evolution 54:563-603.

Thiers, Barbara. 2020. “Index Herbariorum: A Global Directory of Public Herbaria
and Associated Staff.” NYBG Steere Herbarium.
http://sweetgum.nybg.org/science/ih/.
