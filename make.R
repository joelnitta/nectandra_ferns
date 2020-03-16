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

# Load functions and plans
source("code/functions.R")
source("code/plan.R")

# Set cache
nectandra_cache = new_cache("nectandra_cache")

# Specify non-global environment
# to get around captioner modifying global env
# (cf https://github.com/ropensci/drake/issues/749)
envir <- new.env(parent = globalenv())

# Run analyses
make(plan, cache = nectandra_cache, envir = envir) 
