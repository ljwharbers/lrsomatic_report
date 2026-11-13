# Padfoot SV annotation helpers (STUB, WP0 scaffolding — filled in by the SAVANA/Padfoot package).
#
# Plain source()d helpers, no section registration. Sourced by templates/per_sample.qmd and
# tests/testthat/setup.R after R/parse_severus.R and before R/sections/*, so section modules
# (sv.R) may call into it. Use resolve_columns() (R/utils.R) for tolerant column lookup.

suppressPackageStartupMessages({
  library(data.table)
})
