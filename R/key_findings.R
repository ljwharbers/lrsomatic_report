# Clinical overview: key findings + biomarkers (STUB, WP0 scaffolding — body replaced by the
# overview package).
#
# build_overview(variant_table, section_data, all_panels, ascat_pp, tmb_info, counts, ctx)
#   variant_table  FULL small-variant table (build_variant_table()), may be NULL
#   section_data   SECTION_DATA (sv, whatshap, cna, shatterseek, ecdna, reconplot, signatures)
#   all_panels     named list of gene panels (after the YAML round trip: unlist() fields)
#   ascat_pp       list(purity, ploidy) from parse_ascat_purityploidy()
#   tmb_info       compute_tmb() result: list(n_nonsyn, tmb, denominator_mb)
#   counts         list(small_variants = count_small_variants() result, n_sv, n_trans, n_intra_bnd)
#   ctx            the section ctx (see R/sections.R)
# returns
#   list(findings    = data.table(gene, type, alteration, evidence, tier (integer), role,
#                                 panels, vaf, source_id),
#        biomarkers  = list(),
#        per_panel_n = integer())   # named by panel
# Assigned to OVERVIEW in the setup chunk; rendered by templates/sections/_overview.qmd and
# assets/js/key_findings.js. Sourced after R/sections/* by the setup chunk.

suppressPackageStartupMessages({
  library(data.table)
})

empty_findings = function() {
  data.table(gene = character(), type = character(), alteration = character(),
             evidence = character(), tier = integer(), role = character(),
             panels = character(), vaf = numeric(), source_id = character())
}

build_overview = function(variant_table, section_data, all_panels, ascat_pp, tmb_info,
                          counts, ctx) {
  list(findings = empty_findings(), biomarkers = list(), per_panel_n = integer())
}
