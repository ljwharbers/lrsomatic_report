# Copy-number section module: genome-wide CN segments, purity/ploidy, genome-level WGD/LOH
# and per-gene calls for every embedded panel. Parsers live in R/parse_ascat.R.
#
# Contract
#   id "cna"; parse(inputs, section_data, ctx) returns NULL when no CN source, purity/ploidy
#   source or Wakhan solutions table was found, else
#   list(
#     source           = "wakhan" | "ascat" | NA_character_,   # which caller the calls come from
#     purity, ploidy,                                           # of `source` (see precedence below)
#     ploidy_estimated = TRUE when Wakhan segments came without a solutions table and the
#                        ploidy is the length-weighted mean total CN of the autosomes
#     purity_by_source = data.table(source, purity, ploidy, used),   # every source found
#     segments         = data.table(chrom, start, end, cn_a, cn_b, total, minor, conf, sv_ids,
#                                   conf_a, conf_b),
#                        # Wakhan: cn_a/cn_b = HP1/HP2 state, conf = min of the two haplotype
#                        # confidences. ASCAT: cn_a/cn_b = nMajor/nMinor, conf NA. Uncapped.
#     genome           = list(wgd, loh_fraction, wgd_source, major2_fraction, autosomal_bp,
#                             ascat_wgd, ascat_loh, ascat_gi),
#     genes            = data.table(panel, gene, chrom, start, end, min_cn, max_cn, min_minor,
#                                   n_segments, call, partial, loh, source, role, expected),
#                        # call: "amplification" | "homozygous deletion" | "deep deletion" |
#                        # "gain" | "loss" |
#                        # "neutral" | NA (gene span overlaps no segment, n_segments == 0)
#     wakhan_genes     = data.table(chrom, start, end, gene, hp1_state, hp2_state, total, minor),
#     wakhan_genes_dropped = integer(1),   # rows dropped as state > CNA_MAX_PLAUSIBLE_CN
#     panels_without_coords = character(), # panels contributing no coordinate rows
#     solutions        = data.table(),     # wakhan/solutions_ranks.tsv, ordered by rank
#     profiles         = list(list(rank, purity, ploidy, confidence, path, segments)),
#                        # every ranked Wakhan integer profile found (rank 1 first)
#     wakhan_dir       = character(1),     # NA when there is no Wakhan output
#     wakhan_heatmap   = character(1),     # NA unless the heatmap HTML is present
#     ascat_metrics    = list() | NULL     # parse_ascat_metrics()
#   )
#
# Precedence. Segments: Wakhan rank-1 integer_profile.bed, else ASCAT segments_raw.txt.
# Purity/ploidy are those of `source`, so the call thresholds and the segments share one
# ploidy scale; with no segments at all: Wakhan, else ASCAT, else SAVANA. Wakhan's purity is
# `cell_purity` (fraction of tumour cells, the quantity ASCAT's AberrantCellFraction and
# SAVANA's purity estimate), not `dna_purity` (tumour fraction of the DNA).
#
# Calls (COSMIC Cancer Gene Census CNA definition): ploidy <= CNA_PLOIDY_SPLIT -> amplification
# when total CN >= CNA_AMP_CN_LOW, homozygous deletion when total CN == 0; ploidy above it ->
# amplification when CN >= CNA_AMP_CN_HIGH, homozygous deletion when CN < ploidy -
# CNA_LOSS_OFFSET, named "deep deletion" unless the total CN is 0. Otherwise gain/loss against the expected CN: round(ploidy), and for chrX/chrY
# in a male max(1, round(ploidy / 2)) — 1 at any diploid-ish ploidy, 2 after a WGD.
# A gene spanning several segments: amplification on max_cn, homozygous deletion on min_cn,
# then loss, then gain; partial = TRUE when the per-segment calls disagree.

CNA_PLOIDY_SPLIT     = 2.7
CNA_AMP_CN_LOW       = 5
CNA_AMP_CN_HIGH      = 9
CNA_LOSS_OFFSET      = 2.7
CNA_MAX_PLAUSIBLE_CN = 100
CNA_PROFILE_YCAP     = 6      # profile y axis is capped here; segments above get a marker
CNA_PROFILE_MIN_FRAC = 0.001  # profile only: segments narrower than this genome fraction are
                              # drawn widened (centred) so focal events stay visible
CNA_CALL_LEVELS      = c("amplification", "homozygous deletion", "deep deletion", "loss", "gain",
                         "neutral")

# Expected total CN per segment: round(ploidy); halved (min 1) on chrX/chrY in a male
cna_expected_cn = function(chrom, ploidy, sex = "female") {
  male = tolower(trimws(as.character(sex))) %in% c("male", "xy", "m")
  base = round(ploidy)
  if (male) ifelse(chrom %in% c("chrX", "chrY"), pmax(1, round(ploidy / 2)), base)
  else rep(base, length(chrom))
}

# Pure, vectorised: data.table(call, loh). LOH = no minor allele on a non-zero segment where
# two alleles are expected (so a male's hemizygous chrX is not LOH)
cna_call = function(total, minor, ploidy, expected) {
  n = max(length(total), length(minor), length(ploidy), length(expected))
  total = rep_len(as.numeric(total), n); minor = rep_len(as.numeric(minor), n)
  ploidy = rep_len(as.numeric(ploidy), n); expected = rep_len(as.numeric(expected), n)
  low = ploidy <= CNA_PLOIDY_SPLIT
  amp = total >= ifelse(low, CNA_AMP_CN_LOW, CNA_AMP_CN_HIGH)
  hd  = ifelse(low, total == 0, total < ploidy - CNA_LOSS_OFFSET)
  # Above the ploidy split the CGC threshold can fire on CN >= 1: that is a deep, not a homozygous, deletion
  call = fifelse(amp, "amplification",
         fifelse(hd, fifelse(total == 0, "homozygous deletion", "deep deletion"),
         fifelse(total < expected, "loss",
         fifelse(total > expected, "gain", "neutral"))))
  call[is.na(total) | is.na(ploidy) | is.na(expected)] = NA_character_
  loh = !is.na(minor) & !is.na(total) & minor == 0 & total > 0 & !is.na(expected) & expected >= 2
  data.table(call = call, loh = loh)
}

# Genome-level metrics over the autosomes, as fractions of the autosomal length covered by
# segments: WGD when the major allele is >= 2 over more than half of it
cna_genome_metrics = function(segments) {
  na = list(wgd = NA, loh_fraction = NA_real_, major2_fraction = NA_real_, autosomal_bp = 0)
  if (is.null(segments) || nrow(segments) == 0) return(na)
  s = segments[chrom %in% paste0("chr", 1:22) & !is.na(total)]
  if (nrow(s) == 0) return(na)
  len = pmax(0, s$end - s$start + 1)
  covered = sum(len)
  if (covered <= 0) return(na)
  major = pmax(s$cn_a, s$cn_b)
  major2 = sum(len[major >= 2]) / covered
  list(wgd = major2 > 0.5,
       loh_fraction = sum(len[s$minor == 0 & s$total > 0]) / covered,
       major2_fraction = major2, autosomal_bp = covered)
}

.empty_cna_genes = function() {
  data.table(panel = character(), gene = character(), chrom = character(),
             start = numeric(), end = numeric(), min_cn = numeric(), max_cn = numeric(),
             min_minor = numeric(), n_segments = integer(), call = character(),
             partial = logical(), loh = logical(), source = character(),
             role = character(), expected = numeric(), span_from = character())
}

# Collapse per-segment calls to one gene call: most severe present, in CNA_CALL_LEVELS order
.cna_collapse_call = function(calls) {
  calls = unique(calls[!is.na(calls)])
  if (length(calls) == 0) return(NA_character_)
  CNA_CALL_LEVELS[min(match(calls, CNA_CALL_LEVELS))]
}

# Per-gene calls for every panel. Coordinate rows (panel_gene_intervals(p, "all")) are
# overlapped with `segments`. A panel symbol without coordinates takes the gene's span from
# another embedded list that has one (span_from = that list), and only failing that Wakhan's
# own gene file (`wakhan_genes`, symbol match, called at `wakhan_ploidy`; span_from = "Wakhan
# segment": its coordinates are the segment's, not the gene's). span_from is NA for a list's
# own coordinates. Returns list(genes, panels_without_coords).
cna_gene_calls = function(segments, panels, ploidy, sex, source,
                          wakhan_genes = NULL, wakhan_ploidy = ploidy) {
  out = list()
  no_coords = character(0)
  have_segs = !is.null(segments) && nrow(segments) > 0 && !is.na(ploidy)
  if (have_segs) {
    segs = segments[!is.na(total), .(chrom, start, end, total, minor)]
    setkey(segs, chrom, start, end)
  }
  # Every embedded list's gene spans, for symbols another list leaves without coordinates
  spans = rbindlist(lapply(names(panels), function(nm) {
    v = panel_gene_intervals(panels[[nm]], "all")
    if (nrow(v)) v[, from := nm]
  }))
  for (nm in names(panels)) {
    p = panels[[nm]]
    iv = panel_gene_intervals(p, "all")
    if (nrow(iv) == 0) no_coords = c(no_coords, nm)
    iv[, from := NA_character_]
    # Symbols without coordinates here: borrow the span another list gives the gene
    lend = setdiff(unique(as.character(unlist(p$genes))), iv$gene)
    if (length(lend) > 0 && nrow(spans) > 0) {
      b = spans[toupper(gene) %in% toupper(lend) & from != nm]
      if (nrow(b)) {
        b[, gk := toupper(gene)]
        b = b[b[, .(from1 = from[1]), by = gk], on = "gk"][from == from1]
        b[, gene := lend[match(gk, toupper(lend))]]
        iv = rbind(iv, b[, .(chrom, start, end, gene, from)])
      }
    }
    if (nrow(iv) > 0 && have_segs) {
      g = iv[, .(gstart = min(start), gend = max(end), from = from[1]), by = .(gene, chrom)]
      ov = foverlaps(g, segs, by.x = c("chrom", "gstart", "gend"), type = "any",
                     nomatch = NA)
      ov[, expected := cna_expected_cn(chrom, ploidy, sex)]
      cc = cna_call(ov$total, ov$minor, ploidy, ov$expected)
      ov[, `:=`(seg_call = cc$call, seg_loh = cc$loh)]
      agg = ov[, {
        ok = !is.na(total)
        calls = unique(seg_call[ok])
        list(min_cn = if (any(ok)) min(total[ok]) else NA_real_,
             max_cn = if (any(ok)) max(total[ok]) else NA_real_,
             min_minor = if (any(ok & !is.na(minor))) min(minor[ok], na.rm = TRUE) else NA_real_,
             n_segments = sum(ok),
             call = .cna_collapse_call(seg_call[ok]),
             partial = length(calls[!is.na(calls)]) > 1,
             loh = any(seg_loh[ok]),
             expected = expected[1])
      }, by = .(gene, chrom, start = gstart, end = gend, span_from = from)]
      agg[, `:=`(panel = nm, source = source)]
      out[[length(out) + 1L]] = agg
    }
    # Symbols without coordinates: Wakhan's gene file, if it has the symbol
    syms = setdiff(unique(as.character(unlist(p$genes))), iv$gene)
    if (length(syms) > 0 && !is.null(wakhan_genes) && nrow(wakhan_genes) > 0 &&
        !is.na(wakhan_ploidy)) {
      w = wakhan_genes[toupper(gene) %in% toupper(syms)]
      w = w[!duplicated(toupper(gene))]
      if (nrow(w) > 0) {
        w[, expected := cna_expected_cn(chrom, wakhan_ploidy, sex)]
        cc = cna_call(w$total, w$minor, wakhan_ploidy, w$expected)
        out[[length(out) + 1L]] = data.table(
          gene = w$gene, chrom = w$chrom, start = w$start, end = w$end,
          min_cn = w$total, max_cn = w$total, min_minor = w$minor, n_segments = 1L,
          call = cc$call, partial = FALSE, loh = cc$loh, expected = w$expected,
          panel = nm, source = "wakhan-genes", span_from = "Wakhan segment")
      }
    }
  }
  genes = if (length(out) > 0) rbindlist(out, fill = TRUE) else .empty_cna_genes()
  genes[, role := NA_character_]
  for (nm in names(panels)) {
    r = panel_roles(panels[[nm]])
    if (nrow(r) > 0) genes[panel == nm, role := r$role[match(gene, r$gene)]]
  }
  genes[, n_segments := as.integer(n_segments)]
  setcolorder(genes, names(.empty_cna_genes()))
  list(genes = genes[], panels_without_coords = no_coords)
}

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

cna_locate = function(sample_dir, sample_id) {
  hits = list.files(sample_dir,
                    pattern = paste0("\\.segments_raw\\.txt$|\\.purityploidy\\.txt$|",
                                     "\\.metrics\\.txt$|_fitted_purity_ploidy\\.tsv$|",
                                     "^solutions_ranks\\.tsv$"),
                    recursive = TRUE, full.names = TRUE)
  # Prefer a file named for the sample; else the first hit
  pick = function(pattern) {
    h = hits[grepl(pattern, basename(hits))]
    if (length(h) == 0) return(NULL)
    own = h[startsWith(basename(h), paste0(sample_id, ".")) |
            startsWith(basename(h), paste0(sample_id, "_"))]
    if (length(own) > 0) own[1] else h[1]
  }
  seg = pick("\\.segments_raw\\.txt$")
  pp  = pick("\\.purityploidy\\.txt$")
  # metrics.txt is a generic name: only ASCAT's, i.e. beside and prefixed like the
  # purityploidy (or segments) file
  anchor = if (!is.null(pp)) pp else seg
  metrics = if (!is.null(anchor)) {
    f = file.path(dirname(anchor),
                  sub("\\.(purityploidy|segments_raw)\\.txt$", ".metrics.txt", basename(anchor)))
    if (file.exists(f)) f else NULL
  }

  # Wakhan sits at a fixed sample_dir/wakhan (real output and staged layout alike); a
  # solutions table found elsewhere is the fallback
  wdir = file.path(sample_dir, "wakhan")
  if (!dir.exists(wdir)) {
    s = hits[basename(hits) == "solutions_ranks.tsv"]
    wdir = if (length(s) > 0) dirname(s[1]) else NA_character_
  }
  sol = if (!is.na(wdir) && file.exists(file.path(wdir, "solutions_ranks.tsv")))
          file.path(wdir, "solutions_ranks.tsv") else NULL
  heat = if (!is.na(wdir)) Sys.glob(file.path(wdir, "*heatmap_ploidy_purity.html")) else character(0)

  # Ranks to look for: the solutions table's, else rank 1 only
  ranks = data.table(solution_rank = 1L, repository_name = NA_character_)
  if (!is.null(sol)) {
    s = tryCatch(fread(sol, sep = "\t", header = TRUE), error = function(e) NULL)
    if (!is.null(s) && nrow(s) > 0 && "solution_rank" %in% names(s))
      ranks = data.table(solution_rank = as.integer(s$solution_rank),
                         repository_name = if ("repository_name" %in% names(s))
                           as.character(s$repository_name) else NA_character_)
  }
  setorder(ranks, solution_rank)
  profiles = list()
  if (!is.na(wdir)) for (i in seq_len(nrow(ranks))) {
    r = ranks$solution_rank[i]; rn = ranks$repository_name[i]
    bed = locate_wakhan_solution_file(wdir, r, rn, "integer_profile.bed")
    if (is.null(bed)) next
    profiles[[length(profiles) + 1L]] = list(
      rank = r, bed = bed,
      genes = locate_wakhan_solution_file(wdir, r, rn,
                                          file.path("genes", "genes_copynumber_states.bed")))
  }

  list(ascat_segments = seg, ascat_purityploidy = pp, ascat_metrics = metrics,
       savana_purity = pick("_fitted_purity_ploidy\\.tsv$"),
       wakhan = list(dir = wdir, solutions = sol,
                     heatmap = if (length(heat) > 0) heat[1] else NA_character_,
                     profiles = profiles))
}

# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------

cna_parse = function(inputs, section_data, ctx) {
  sex = if (!is.null(ctx$sex)) as.character(unlist(ctx$sex))[1] else "female"
  w = inputs$wakhan

  solutions = parse_wakhan_solutions(w$solutions)
  sol_row = function(rank) {
    if (is.null(solutions)) return(NULL)
    r = solutions[solution_rank == rank]
    if (nrow(r) == 0) NULL else r[1]
  }
  profiles = Filter(Negate(is.null), lapply(w$profiles, function(p) {
    seg = parse_wakhan_integer_profile(p$bed)
    if (is.null(seg) || nrow(seg) == 0) return(NULL)
    r = sol_row(p$rank)
    list(rank = p$rank,
         purity = if (!is.null(r) && "cell_purity" %in% names(r)) as.numeric(r$cell_purity) else NA_real_,
         ploidy = if (!is.null(r)) as.numeric(r$ploidy) else NA_real_,
         confidence = if (!is.null(r) && "confidence" %in% names(r)) as.numeric(r$confidence) else NA_real_,
         path = p$bed, genes_path = p$genes, segments = seg)
  }))
  top = if (length(profiles) > 0) profiles[[1]] else NULL

  ascat_seg = parse_ascat_integer_segments(inputs$ascat_segments)
  ascat_pp  = if (!is.null(inputs$ascat_purityploidy))
                parse_ascat_purityploidy(inputs$ascat_purityploidy) else NULL
  ascat_metrics = parse_ascat_metrics(inputs$ascat_metrics)
  savana = parse_savana_purity(inputs$savana_purity)

  # Purity/ploidy of every source found
  wk_pp = if (!is.null(solutions)) {
    r = solutions[1]
    list(purity = if ("cell_purity" %in% names(r)) as.numeric(r$cell_purity) else NA_real_,
         ploidy = as.numeric(r$ploidy))
  }
  pbs = rbindlist(list(
    if (!is.null(wk_pp)) data.table(source = "wakhan", purity = wk_pp$purity, ploidy = wk_pp$ploidy),
    if (!is.null(ascat_pp) && !(is.na(ascat_pp$purity) && is.na(ascat_pp$ploidy)))
      data.table(source = "ascat", purity = ascat_pp$purity, ploidy = ascat_pp$ploidy),
    if (!is.null(savana)) data.table(source = "savana", purity = savana$purity, ploidy = savana$ploidy)
  ))
  if (nrow(pbs) == 0) pbs = data.table(source = character(), purity = numeric(), ploidy = numeric())

  source = if (!is.null(top)) "wakhan" else
           if (!is.null(ascat_seg) && nrow(ascat_seg) > 0) "ascat" else NA_character_
  segments = if (identical(source, "wakhan")) top$segments else
             if (identical(source, "ascat")) ascat_seg else NULL

  if (is.na(source) && nrow(pbs) == 0 && is.null(solutions)) return(NULL)

  # Purity/ploidy of `source`; with no segments, the first source in precedence order
  pp_src = if (!is.na(source)) source else if (nrow(pbs) > 0) pbs$source[1] else NA_character_
  purity = if (!is.na(pp_src) && pp_src %in% pbs$source) pbs[source == pp_src]$purity[1] else NA_real_
  ploidy = if (!is.na(pp_src) && pp_src %in% pbs$source) pbs[source == pp_src]$ploidy[1] else NA_real_
  if (identical(source, "wakhan") && !is.null(top) && !is.na(top$ploidy)) {
    purity = top$purity; ploidy = top$ploidy
  }
  ploidy_estimated = FALSE
  if (!is.na(source) && is.na(ploidy)) {
    a = segments[chrom %in% paste0("chr", 1:22) & !is.na(total)]
    len = pmax(0, a$end - a$start + 1)
    if (sum(len) > 0) { ploidy = round(sum(a$total * len) / sum(len), 2); ploidy_estimated = TRUE }
  }
  pbs[, used := source == pp_src]

  genome = cna_genome_metrics(segments)
  genome$wgd_source = if (!is.na(source)) source else NA_character_
  genome$ascat_wgd = if (!is.null(ascat_metrics)) ascat_metrics$wgd else NA_real_
  genome$ascat_loh = if (!is.null(ascat_metrics)) ascat_metrics$loh else NA_real_
  genome$ascat_gi  = if (!is.null(ascat_metrics)) ascat_metrics$gi else NA_real_

  wg = if (!is.null(top)) parse_wakhan_genes(top$genes_path, CNA_MAX_PLAUSIBLE_CN) else NULL
  wakhan_genes = if (!is.null(wg)) wg$genes else NULL
  wk_ploidy = if (!is.null(top) && !is.na(top$ploidy)) top$ploidy else ploidy

  panels = if (!is.null(ctx$all_panels)) ctx$all_panels else list()
  gc = cna_gene_calls(segments, panels, ploidy, sex, source, wakhan_genes, wk_ploidy)

  list(
    source = source, purity = purity, ploidy = ploidy, ploidy_estimated = ploidy_estimated,
    purity_by_source = pbs,
    segments = if (!is.null(segments)) segments else .empty_cn_segments(),
    genome = genome,
    genes = gc$genes,
    wakhan_genes = if (!is.null(wakhan_genes)) wakhan_genes else data.table(),
    wakhan_genes_dropped = if (!is.null(wg)) as.integer(wg$n_dropped) else 0L,
    panels_without_coords = gc$panels_without_coords,
    solutions = if (!is.null(solutions)) solutions else data.table(),
    profiles = profiles,
    wakhan_dir = if (!is.null(w$dir)) w$dir else NA_character_,
    wakhan_heatmap = if (!is.null(w$heatmap)) w$heatmap else NA_character_,
    ascat_metrics = ascat_metrics
  )
}

# ---------------------------------------------------------------------------
# Static haplotype-resolved profile (replaces the embedded Wakhan Plotly HTML)
# ---------------------------------------------------------------------------

# ggplot of cn_a/cn_b along a concatenated genome axis; NULL when nothing is plottable.
# `labels` names the two tracks (HP1/HP2 for Wakhan, major/minor for ASCAT).
cna_profile_plot = function(segments, chromosomes, chrom_lens, labels = c("HP1", "HP2"),
                            cap = CNA_PROFILE_YCAP) {
  if (is.null(segments) || nrow(segments) == 0) return(NULL)
  chroms = chromosomes[chromosomes %in% names(chrom_lens)]
  if (length(chroms) == 0) return(NULL)
  lens = as.numeric(chrom_lens[chroms])
  off  = setNames(c(0, cumsum(lens)[-length(lens)]), chroms)
  s = segments[chrom %in% chroms]
  if (nrow(s) == 0) return(NULL)

  dodge = 0.09
  long = rbindlist(list(
    s[, .(chrom, start, end, cn = cn_a, conf = conf_a, track = labels[1], dy = dodge)],
    s[, .(chrom, start, end, cn = cn_b, conf = conf_b, track = labels[2], dy = -dodge)]
  ))[!is.na(cn)]
  long[, `:=`(x = start + off[chrom], xend = end + off[chrom])]
  # Widen focal segments to a visible minimum and draw them last (on top of their neighbours)
  min_w = CNA_PROFILE_MIN_FRAC * sum(lens)
  long[, w := xend - x]
  long[w < min_w, `:=`(x = (x + xend) / 2 - min_w / 2, xend = (x + xend) / 2 + min_w / 2)]
  long[, `:=`(x = pmax(x, 0), xend = pmin(xend, sum(lens)))]
  setorder(long, -w)
  long[, y := pmin(cn, cap) + dy]
  long[is.na(conf), conf := 1]
  above = long[cn > cap]
  above[, xm := (x + xend) / 2]

  bands = data.table(xmin = off, xmax = off + lens, i = seq_along(chroms))[i %% 2 == 0]
  mids  = off + lens / 2
  cols  = setNames(c("#b3402f", "#0d5c75"), labels)

  p = ggplot2::ggplot() +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "#f1eee6", colour = NA) +
    ggplot2::geom_segment(data = long, ggplot2::aes(x = x, xend = xend, y = y, yend = y,
                                                    colour = track, alpha = conf),
                          linewidth = 1.3, lineend = "butt") +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_alpha_continuous(range = c(0.25, 1), limits = c(0, 1), breaks = c(0.25, 0.5, 1),
                                    name = "confidence") +
    ggplot2::guides(alpha = ggplot2::guide_legend(override.aes = list(colour = "#5b6770"))) +
    ggplot2::scale_x_continuous(breaks = mids, labels = sub("^chr", "", chroms),
                                expand = c(0, 0), limits = c(0, sum(lens))) +
    ggplot2::scale_y_continuous(breaks = 0:cap, labels = c(0:(cap - 1), paste0(cap, "+")),
                                limits = c(-0.35, cap + 0.85), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Copy number") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "top", legend.justification = "right",
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(2, "pt"),
                   legend.title = ggplot2::element_text(size = 8),
                   legend.text = ggplot2::element_text(size = 8),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 7),
                   panel.grid.major.y = ggplot2::element_line(colour = "#e6e2d8", linewidth = 0.25))
  if (nrow(above) > 0) {
    p = p +
      ggplot2::geom_point(data = above, ggplot2::aes(x = xm, y = cap + 0.38, colour = track),
                          shape = 17, size = 1.6, show.legend = FALSE) +
      ggplot2::geom_text(data = above, ggplot2::aes(x = xm, y = cap + 0.68, label = cn, colour = track),
                         size = 2.3, show.legend = FALSE)
  }
  p
}

register_section(list(
  id     = "cna",
  title  = "Copy number",
  locate = cna_locate,
  parse  = cna_parse
))
