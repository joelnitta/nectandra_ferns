# Install packages to a docker image with renv

# Updates to packages can be made as follows:
# 1. Volume-mount this repo to the image.
# 2. Set repo locations: 
# ```
# my_repos <- BiocManager::repositories()
# my_repos["CRAN"] <- "https://cran.rstudio.com/"
# options(repos = my_repos)
# ```
# 2. Install new packages to `/renv`.
# 3. Run `renv::snapshot(type = "simple", library = "/renv")`.

### Initialize renv ###

# Install renv
install.packages("remotes", repos = "https://cran.rstudio.com/")
remotes::install_github("rstudio/renv")

# Initialize renv, but don't let it try to find packages to install itself.
renv::init(
  bare = TRUE,
  force = TRUE,
  restart = FALSE)

renv::activate()

### Setup repositories ###

# Install packages that install packages.
install.packages("remotes", repos = "https://cran.rstudio.com/")
install.packages("BiocManager", repos = "https://cran.rstudio.com/")

# Specify repositories so they get included in
# renv.lock file.
my_repos <- BiocManager::repositories()
my_repos["CRAN"] <- "https://cran.rstudio.com/"
options(repos = my_repos)

### Install packages ###

# All packages will be installed to
# the project-specific library.

# Install CRAN packages
cran_packages <- c(
  "assertr",
  "broom",
  "conflicted",
  "cowplot",
  "drake",
  "ggrepel",
  "ggridges",
  "here",
  "iNEXT",
  "ips",
  "janitor",
  "kableExtra",
  "magick",
  "pryr",
  "raster",
  "readxl",
  "writexl",
  "rmarkdown",
  "scico",
  "scales",
  "tidyverse")

install.packages(cran_packages)

### Install bioconductor packages ###
bioc_packages <- c("ggtree")
BiocManager::install(bioc_packages)

# Install github packages
github_packages <- c(
  "joelnitta/jntools",
  "joelnitta/gbfetch",
  "thomasp85/patchwork",
  "rstudio/gt",
  "adletaw/captioner",
  "joelnitta/taxastand")

remotes::install_github(github_packages)

### Take snapshot ###

renv::snapshot(type = "simple")
