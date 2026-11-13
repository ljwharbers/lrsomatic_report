# Copy-number section module (STUB, WP0 scaffolding — body replaced by the CNA package).
#
# Contract
#   id "cna"; parse(inputs, section_data, ctx) returns NULL when nothing is available, else
#   list(
#     source           = "wakhan" | "ascat" | NA_character_,   # which caller the calls come from
#     purity, ploidy,                                           # of `source`
#     purity_by_source = data.table(source, purity, ploidy),
#     segments         = data.table(chrom, start, end, cn_a, cn_b, total, minor, conf, sv_ids),
#     genome           = list(wgd, loh_fraction, wgd_source),
#     genes            = data.table(panel, gene, chrom, start, end, min_cn, max_cn, min_minor,
#                                   n_segments, call, partial, loh, source, role),
#                        # call: "amplification" | "homozygous deletion" | "gain" | "loss" | "neutral"
#     wakhan_genes     = data.table(),
#     solutions        = data.table(),
#     wakhan_dir       = character(1)
#   )
#   Gene intervals come from panel_gene_intervals(panel, "all") (R/utils.R), roles from
#   panel_roles(panel). Rendered by templates/sections/_cna_genes.qmd (end of _ascat.qmd).

register_section(list(
  id     = "cna",
  title  = "Copy number",
  locate = function(sample_dir, sample_id) list(),
  parse  = function(inputs, section_data, ctx) NULL
))
