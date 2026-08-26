# Assumes tests run from repo root (default for testthat::test_dir())
library(testthat)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Source R modules from repo root
# Use absolute paths to handle testthat::test_dir sourcing context
repo_root = dirname(dirname(getwd()))
source(file.path(repo_root, "R/utils.R"))
source(file.path(repo_root, "R/references.R"))
source(file.path(repo_root, "R/parse_smallvariants.R"))
source(file.path(repo_root, "R/parse_severus.R"))
# The breakend circos is drawn client-side now, so R/circos_bnd.R only selects and
# serialises data — no circlize, nothing to skip. (R/circos.R still needs circlize for
# the genome-wide plot, but nothing under test sources it.)
source(file.path(repo_root, "R/sections.R"))
source(file.path(repo_root, "R/circos_bnd.R"))
