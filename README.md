# nectandra_ferns

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "A taxonomic and molecular survey of the pteridophytes of the Nectandra Reserve, Costa Rica".

All code is in [R](https://cran.r-project.org/). The [drake package](https://ropensci.github.io/drake/) is used to manage the workflow. To run all analyses and generate the manuscript, [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) and run `make.R`.

## NCBI Entrez API key

This analysis queries GenBank and downloads DNA sequences. It is recommended to register a free [NCBI
Entrez API
key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) and set this as the `ENTREZ_KEY` environmental variable before running analyses, or warnings will be issued. See [`taxize::use_entrez()`](https://www.rdocumentation.org/packages/taxize/versions/0.9.94/topics/key_helpers) and [``taxize::`taxize-authentication`()``](https://www.rdocumentation.org/packages/taxize/versions/0.9.94/topics/taxize-authentication) for more info on how to set this up.

## Reproducible analysis with Docker

`make.R` requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/nectandra) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/) and clone this repository.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

Enter the container:

```
docker exec -it nectandra_ferns_analysis_1 bash
```

Inside the container, run `make.R`:

```
Rscript make.R
```

You will see the targets being built by `drake`, and the final manuscript should be compiled at the end as `nectandra_ferns.pdf` and `nectandra_ferns.docx` in the `results` folder. Other figure and table files will also be compiled.

When it's finished, exit the container and take it down:

```
exit
docker-compose down
```
