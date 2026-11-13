# ReConPlot figures section module (STUB, WP0 scaffolding).
#
# Contract
#   id "reconplot"; parse(inputs, section_data, ctx) returns NULL or
#   list(combos = list(<combo> = list(genome = <path>,
#                                     chroms = c(chr1 = <path>, ...))))   # named chr vector
#   Rendered by templates/sections/_reconplot.qmd (included at the end of _sv.qmd).

register_section(list(
  id     = "reconplot",
  title  = "ReConPlot",
  locate = function(sample_dir, sample_id) list(),
  parse  = function(inputs, section_data, ctx) NULL
))
