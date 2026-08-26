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
# circlize is a hard dependency of the breakend circos; the tests that need it skip
# themselves, so a machine without it still runs the rest of the suite.
if (requireNamespace("circlize", quietly = TRUE)) {
  source(file.path(repo_root, "R/circos_bnd.R"))
}
