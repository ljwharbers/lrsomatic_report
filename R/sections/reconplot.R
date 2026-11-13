# ReConPlot figures section module.
#
# Inputs: reconplot/<combo>/genome_wide/<id>_genome_wide.png and
#         reconplot/<combo>/per_chromosome/<id>_chr*.png, combo in severus_ascat, severus_wakhan,
#         savana (any other combo directory is kept, after those). PDFs, focus/ figures and the
#         harmonised TSVs are ignored.
#
# Contract
#   id "reconplot"; parse(inputs, section_data, ctx) returns NULL (no figure at all) or
#   list(combos = list(<combo> = list(genome = <path or NULL>,
#                                     chroms = c(chr1 = <path>, ...))))   # named chr vector
#   chroms is in display order (ctx$chromosomes, unknown contigs last). Rendered by
#   templates/sections/_reconplot.qmd (included at the end of _sv.qmd).

RECONPLOT_COMBOS = c("severus_ascat", "severus_wakhan", "savana")
RECONPLOT_COMBO_LABELS = c(severus_ascat = "Severus + ASCAT", severus_wakhan = "Severus + Wakhan",
                           savana = "SAVANA")

register_section(list(
  id     = "reconplot",
  title  = "ReConPlot",

  locate = function(sample_dir, sample_id) {
    # Recursive file listing sees the staged symlinked files (not symlinked directories)
    pngs = list.files(sample_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
    m = regmatches(pngs, regexec("(^|/)reconplot/([^/]+)/(genome_wide|per_chromosome)/([^/]+)\\.png$",
                                 pngs))
    ok = lengths(m) > 0
    if (!any(ok)) return(list(files = data.table(combo = character(), kind = character(),
                                                 file = character(), path = character()),
                              sample_id = sample_id))
    parts = do.call(rbind, m[ok])
    list(files = data.table(combo = parts[, 3], kind = parts[, 4], file = parts[, 5],
                            path = pngs[ok]),
         sample_id = sample_id)
  },

  parse  = function(inputs, section_data, ctx) {
    f = inputs$files
    if (is.null(f) || nrow(f) == 0) return(NULL)
    sid = inputs$sample_id

    # Several samples' figures in one tree: keep this sample's, where it can be told apart
    own = function(d, suffix) {
      hit = d[file == paste0(sid, suffix)]
      if (nrow(hit) > 0) hit else d
    }

    combos = unique(f$combo)
    combos = c(intersect(RECONPLOT_COMBOS, combos), sort(setdiff(combos, RECONPLOT_COMBOS)))
    chrom_order = if (!is.null(ctx$chromosomes)) as.character(unlist(ctx$chromosomes)) else character(0)

    out = list()
    for (cb in combos) {
      g = own(f[combo == cb & kind == "genome_wide" & grepl("_genome_wide$", file)],
              "_genome_wide")
      genome = if (nrow(g) > 0) g$path[order(nchar(g$path))][1] else NULL

      pc = f[combo == cb & kind == "per_chromosome" & grepl("_chr[^_]+$", file)]
      if (nrow(pc) > 0) {
        pc[, chrom := sub("^.*_(chr[^_]+)$", "\\1", file)]
        mine = pc[startsWith(file, paste0(sid, "_chr"))]
        if (nrow(mine) > 0) pc = mine
        pc = pc[!duplicated(chrom)]
        rank = match(pc$chrom, chrom_order)
        # Off-list contigs after the known ones, in natural order (chr2 before chr10)
        nat = suppressWarnings(as.numeric(sub("^chr", "", pc$chrom)))
        pc = pc[order(is.na(rank), rank, is.na(nat), nat, chrom)]
        chroms = setNames(pc$path, pc$chrom)
      } else {
        chroms = setNames(character(0), character(0))
      }

      if (!is.null(genome) || length(chroms) > 0)
        out[[cb]] = list(genome = genome, chroms = chroms)
    }
    if (length(out) == 0) return(NULL)
    list(combos = out)
  }
))
