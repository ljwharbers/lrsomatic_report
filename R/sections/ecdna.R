# ecDNA / amplicon section module (STUB, WP0 scaffolding).
#
# Contract
#   id "ecdna"; parse(inputs, section_data, ctx) returns NULL (stub) or
#   list(
#     status    = "not_run" | "none_detected" | "detected",
#     amplicons = data.table(amplicon, class, ecdna, bfb, location, oncogenes, max_cn),
#     genes     = data.table(amplicon, feature, gene, gene_cn),
#     plots     = character()   # image paths
#   )
#   Rendered by templates/sections/_instability.qmd.

register_section(list(
  id     = "ecdna",
  title  = "ecDNA",
  locate = function(sample_dir, sample_id) list(),
  parse  = function(inputs, section_data, ctx) NULL
))
