# Section-module contract: id, title, locate(sample_dir, sample_id), parse(inputs, section_data[, ctx]); see CLAUDE.md
#
# parse() may take an optional third argument `ctx`. The setup chunk of
# templates/per_sample.qmd passes it only when length(formals(s$parse)) >= 3, so a
# two-argument parse() (sv.R, whatshap.R) keeps working unchanged. `ctx` is a named list:
#   sample_id, sample_dir  the render's identifiers (same values locate() received)
#   reference              "t2t" | "hg38"
#   sex                    "male" | "female"
#   all_panels             named list of gene panels (load_gene_panel() objects, after the
#                          Quarto execute_params YAML round trip — unlist() vector fields)
#   chromosomes            plotted chromosomes for this sex, in display order
#   chrom_lens             named numeric, chromosome -> length (reference-specific)
#   cytobands              data.table from load_cytobands()
#   embed                  "slim" | "full" (--embed)
#   repo_dir               repository root, for assets/
#   outputs                locate_outputs() list; prefer the section's own locate() for
#                          anything it owns, this is for shared inputs (e.g. ASCAT files)
#
# Registration order = parse order, and later sections read earlier ones through
# `section_data`. Files under R/sections/ are therefore sourced in SECTION_ORDER, not
# alphabetically; a file not named there is sourced after them, alphabetically.

SECTIONS = list()

SECTION_ORDER = c("sv", "whatshap", "cna", "shatterseek", "ecdna", "reconplot", "signatures")

register_section = function(descriptor) {
  SECTIONS[[descriptor$id]] <<- descriptor
}

# Section module files in parse order: SECTION_ORDER first, any others alphabetically after
section_files = function(dir) {
  files = sort(list.files(dir, pattern = "\\.R$", full.names = TRUE))
  ids   = sub("\\.R$", "", basename(files))
  rank  = match(ids, SECTION_ORDER)
  rank[is.na(rank)] = length(SECTION_ORDER) + 1L
  files[order(rank, ids)]
}

# Run one section's parse(), passing ctx only to a parse() that declares a third argument
run_section_parse = function(s, inputs, section_data, ctx) {
  if (length(formals(s$parse)) >= 3) s$parse(inputs, section_data, ctx)
  else s$parse(inputs, section_data)
}

# Quiet "nothing to show" notice; `warn` marks a genuine failure
section_notice = function(msg, warn = FALSE) {
  tags$div(class = if (warn) "section-notice section-notice--warn" else "section-notice", msg)
}

# Collapsed provenance/caveats under a table; `summary` names the contents
table_details = function(..., summary = "Details") {
  tags$details(
    class = "table-details",
    tags$summary(summary),
    tags$div(class = "table-footnote", ...)
  )
}
