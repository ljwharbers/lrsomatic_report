# Mutational signatures section module (STUB, WP0 scaffolding).
#
# Contract
#   id "signatures"; parse(inputs, section_data, ctx) returns NULL or
#   list(
#     cosmic_version,
#     contexts = list(
#       SBS = list(activities = data.table(signature, mutations, pct, aetiology, artefact),
#                  stats      = list(total, cosine, l1_pct, l2_pct, kl, correlation),
#                  matrix     = data.table(channel, group, count),
#                  flags      = character()),
#       DBS = <same shape>,
#       ID  = <same shape>),
#     dominant = data.table(context, signature, pct, aetiology, artefact)
#   )
#   Rendered by templates/sections/_signatures.qmd.

register_section(list(
  id     = "signatures",
  title  = "Mutational signatures",
  locate = function(sample_dir, sample_id) list(),
  parse  = function(inputs, section_data, ctx) NULL
))
