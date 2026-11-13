# ShatterSeek chromothripsis section module (STUB, WP0 scaffolding).
#
# Contract
#   id "shatterseek"; parse(inputs, section_data, ctx) returns NULL or
#   list(
#     status = "not_run" | "ok",
#     combos = named list of data.table,            # one per SV-caller/CN-source combination
#     high   = data.table(combo, chrom, confidence)  # high-confidence chromothripsis calls
#   )
#   Rendered by templates/sections/_instability.qmd.

register_section(list(
  id     = "shatterseek",
  title  = "Chromothripsis (ShatterSeek)",
  locate = function(sample_dir, sample_id) list(),
  parse  = function(inputs, section_data, ctx) NULL
))
