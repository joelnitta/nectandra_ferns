# make.R

# Master script for running analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Set random number seed for reproducibility
set.seed(9130)

# Set working directory
setwd(here::here())

# Load packages
source("code/packages.R")

# Update drake settings
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Load functions
source("code/functions.R")
# Load plan
source("code/plan.R")
# Load captions
source("ms/captions.R")

# Set cache
nectandra_cache = new_cache("nectandra_cache")
options(rstudio_drake_cache = nectandra_cache)

# Run analyses
make(plan, cache = nectandra_cache)
