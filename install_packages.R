# Install packages to a docker image with renv

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

# Specify repositories so they get included in
# renv.lock file.
my_repos <- c()
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
  "drake",
  "ggrepel",
  "ggridges",
  "here",
  "iNEXT",
  "ips",
  "janitor",
  "kableExtra",
  "raster",
  "readxl",
  "rmarkdown",
  "scico",
  "tidyverse")

install.packages(cran_packages)

# Install github packages
github_packages <- c(
  "joelnitta/jntools",
  "thomasp85/patchwork",
  "rstudio/gt",
  "r-lib/scales")

remotes::install_github(github_packages)

### Take snapshot ###

renv::snapshot(type = "simple")
