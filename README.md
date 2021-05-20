# nectandra_ferns

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "A taxonomic and molecular survey of the pteridophytes of the Nectandra Reserve, Costa Rica". PLOS ONE 2020 https://doi.org/10.1371/journal.pone.0241231

All code is in [R](https://cran.r-project.org/). The [drake package](https://ropensci.github.io/drake/) is used to manage the workflow. 

To run all analyses and generate the manuscript:

1. [Clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
2. If you haven't already, [set up an NCBI Entrez API key](#ncbi-entrez-api-key)
3. [Download the data](#data)
4. [Run `make.R` in the provided docker container](#reproducible-analysis-with-docker)

## NCBI Entrez API key

This analysis queries GenBank and downloads DNA sequences. It is recommended to register a free [NCBI
Entrez API
key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) and set this as the `ENTREZ_KEY` environmental variable before running analyses, or warnings will be issued when the sequences get downloaded. See [`taxize::use_entrez()`](https://www.rdocumentation.org/packages/taxize/versions/0.9.94/topics/key_helpers) and [``taxize::`taxize-authentication`()``](https://www.rdocumentation.org/packages/taxize/versions/0.9.94/topics/taxize-authentication) for more info on how to set this up.

(short version of the `taxize` documentation: make a file called `.Renviron` in the root of this repo, and put this line in it: `ENTREZ_KEY=somekey`, where `somekey` is the API key you set up on [NCBI](https://www.ncbi.nlm.nih.gov/account/))

## Data

There are two small data files included in this repo. The rest of the data are hosted on [Dryad](https://datadryad.org). For each of the links below, click on "Download Dataset", then place the zipped data file (it will have a similar name to the DOI) in the `data` folder of the this repo. This must be done **before** running `make.R`. You can manually unzip the data archives if you want to see the contents, but the code needs the original zipped file in `data/` to run.

- https://doi.org/10.5061/dryad.bnzs7h477
- https://doi.org/10.5061/dryad.4362p32
- https://doi.org/10.5061/dryad.df59g

## Reproducible analysis with Docker

`make.R` requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/nectandra_ferns) to run the code reproducibly. You can [install docker](https://docs.docker.com/install/) from here.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

Run `make.R` inside the container:

```
docker exec nectandra_ferns_analysis_1 Rscript make.R
```

You will see the targets being built by `drake`, and the final manuscript should be compiled at the end as `nectandra_ferns.pdf` and `nectandra_ferns.docx` in the `results/ms` folder. Other figure and table files will also be compiled.

When it's finished, take down the container:

```
docker-compose down
```

## Licenses

- All code in this repository is licensed under the [MIT license](LICENSE)
- [The data](https://doi.org/10.5061/dryad.fqz612jps) are licensed under the [CC0 1.0 Universal Public Domain Dedication license](https://creativecommons.org/publicdomain/zero/1.0/)
- [The paper](https://doi.org/10.1371/journal.pone.0241231) is licensed under the [Creative Commons Attribution License](https://creativecommons.org/licenses/by/4.0/)
- The [Roboto font](https://github.com/google/roboto/) is licensed under the [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)
