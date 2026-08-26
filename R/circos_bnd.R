suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
})

# Breakend circos: a second, BND-only plot drawn over just the chromosomes a breakend
# touches, cross-linked to the SV table so selecting a row highlights its arc.
#
# Unlike draw_circos() in R/circos.R, whose output is embedded as a base64 <img> (a
# replaced element, whose internals no page script can reach), this one is inlined into
# the document and every object in it carries an id. svglite emits no ids or classes of
# its own, so each object is *drawn* in a colour unique to it and the rendered SVG is
# then rewritten — see .tag_bnd_svg().

# Above this the inline SVG stops being worth its bytes: real samples on this cluster
# carry 50-115 arcs after mate collapse, but an unfiltered run can carry ~13k, at
# ~1.4 kB of <polyline> each.
BND_CIRCOS_MAX_LINKS = 2000L

# The BND classes with two loci to draw an arc between; svclass is set by
# parse_severus_somatic_records() / build_sv_table(). A "single breakend" has no partner
# locus and so cannot be drawn at all.
BND_CIRCOS_CLASSES = c("translocation", "intra-chr breakend")

# Sentinel colour families. Colour is the key rather than document order because
# connector lines and BND arcs are both <polyline>: an order-based scheme would have to
# tell them apart by position, and would break the first time circlize reordered
# anything. Each family code is followed by a zero-padded decimal index, so a sentinel
# never contains a hex letter and can be matched without worrying about case.
BND_SENTINEL_FAMILIES = c(link = "10", chrom = "20", gene = "30", line = "40", body = "50")
BND_SENTINEL_RE = "#[1-5]0[0-9]{4}"

.bnd_sentinel = function(family, i) {
  sprintf("#%s%04d", BND_SENTINEL_FAMILIES[[family]], as.integer(i))
}

.attr = function(x) htmltools::htmlEscape(as.character(x), attribute = TRUE)

# ---- Data selection ------------------------------------------------------

# The rearrangements the circos can draw: BNDs with both loci known, on chromosomes the
# report is plotting. Callers pass the result to both bnd_panel_genes() and
# draw_bnd_circos(), so the gene track cannot annotate an arc that was never drawn.
bnd_links = function(sv_table, chromosomes) {
  cols = c("id", "svclass", "chrom_a", "pos_a", "chrom_b", "pos_b")
  empty = data.table(id = character(), svclass = character(),
                     chrom_a = character(), pos_a = integer(),
                     chrom_b = character(), pos_b = integer())
  if (is.null(sv_table) || nrow(sv_table) == 0 || !all(cols %in% names(sv_table)))
    return(empty)

  b = as.data.table(sv_table)[, ..cols]
  b = b[svclass %in% BND_CIRCOS_CLASSES &
        !is.na(chrom_a) & !is.na(pos_a) & !is.na(chrom_b) & !is.na(pos_b) &
        chrom_a %in% chromosomes & chrom_b %in% chromosomes]
  if (nrow(b) == 0) return(empty)

  # Deterministic order, so the sentinel index a row gets is reproducible across renders.
  setorder(b, chrom_a, pos_a, chrom_b, pos_b, id)
  b[]
}

# The panel genes to label: for every coordinate-carrying panel, those whose interval
# falls within `window` of either breakend of a drawn arc.
#
# This is the same test, with the same window, that sv_panel_hits() applies to fill the
# table's panel_hit column — deliberately, because the report only shows a gene here once
# it appears in that column. Reimplementing the predicate would let the two drift.
#
# A symbol-only panel (and the custom paste-in panel) carries no coordinates, so it can
# produce panel_hit labels but no gene bodies.
bnd_panel_genes = function(links, all_panels, window = SV_PANEL_WINDOW_BND) {
  empty = data.table(chrom = character(), start = integer(), end = integer(),
                     gene = character(), panels = character())
  if (is.null(links) || nrow(links) == 0) return(empty)
  if (is.null(all_panels) || length(all_panels) == 0) return(empty)

  q = rbindlist(list(
    data.table(chrom = links$chrom_a, pos = as.numeric(links$pos_a)),
    data.table(chrom = links$chrom_b, pos = as.numeric(links$pos_b))
  ))
  q = unique(q[!is.na(chrom) & !is.na(pos)])
  if (nrow(q) == 0) return(empty)
  q[, `:=`(start = pmax(pos - window, 0), end = pos + window)]

  out = list()
  for (nm in names(all_panels)) {
    iv = panel_intervals(all_panels[[nm]])
    if (is.null(iv) || nrow(iv) == 0) next
    iv = copy(iv)[, `:=`(start = as.numeric(start), end = as.numeric(end))]
    setkey(iv, chrom, start, end)
    ov = data.table::foverlaps(q[, .(chrom, start, end)], iv,
                               by.x = c("chrom", "start", "end"),
                               type = "any", nomatch = NULL)
    if (nrow(ov) == 0) next
    out[[nm]] = unique(ov[, .(chrom, start, end, gene, panel = nm)])
  }
  if (length(out) == 0) return(empty)

  g = rbindlist(out)
  g = g[, .(panels = paste(sort(unique(panel)), collapse = ",")),
        by = .(chrom, start, end, gene)]
  setorder(g, chrom, start, gene)
  g[, `:=`(start = as.integer(start), end = as.integer(end))]
  g[, .(chrom, start, end, gene, panels)]
}

# ---- SVG post-processing -------------------------------------------------

# Swap every sentinel colour for real styling and attach class/data-* attributes.
#
# Works line by line — svglite writes one element per line — and looks each colour up in
# a registry rather than trusting the pattern, so a cytoband grey that merely looks like
# a sentinel (#101010) is left alone. Shapes lose their inline `style` entirely so that
# CSS owns stroke, width and opacity; <text> keeps its font declarations and only has its
# fill swapped for currentColor, which is what makes the labels follow the report's
# dark-mode text colour.
#
# `registry` maps sentinel colour -> list(family, attrs, expect_text). When expect_text is
# given, the element's rendered text must match it or the element is left untagged: that
# is the invariant check on circos.genomicLabels() keeping each label paired with the
# colour it was handed.
.tag_bnd_svg = function(txt, registry) {
  lines = strsplit(txt, "\n", fixed = TRUE)[[1]]
  seen = character(0)
  mismatched = character(0)

  for (i in grep(BND_SENTINEL_RE, lines)) {
    ln = lines[i]
    key = regmatches(ln, regexpr(BND_SENTINEL_RE, ln))
    meta = registry[[key]]
    if (is.null(meta)) next

    if (!is.null(meta$expect_text)) {
      shown = sub("^.*>([^<]*)</text>.*$", "\\1", ln)
      if (!identical(shown, meta$expect_text)) {
        mismatched = c(mismatched, key)
        next
      }
    }

    if (grepl("<text", ln, fixed = TRUE)) {
      ln = sub(key, "currentColor", ln, fixed = TRUE)
      ln = sub("<text", paste0("<text ", meta$attrs), ln, fixed = TRUE)
    } else {
      ln = sub(" style='[^']*'", "", ln)
      ln = sub("<(polyline|polygon|rect|circle|path)",
               paste0("<\\1 ", meta$attrs), ln)
    }

    lines[i] = ln
    seen = c(seen, key)
  }

  # The <style> block svglite puts in <defs> is global page CSS once the SVG is inlined,
  # and its `.svglite polyline { stroke: #000000 }` rule would leak into the document.
  out = paste(lines, collapse = "\n")
  out = gsub(".svglite", ".bnd-svglite", out, fixed = TRUE)
  out = sub("<g class='svglite'>", "<g class='bnd-svglite'>", out, fixed = TRUE)
  out = substring(out, regexpr("<svg", out, fixed = TRUE))

  list(svg = out, seen = seen, mismatched = unique(mismatched))
}

# ---- The plot ------------------------------------------------------------

# Draw the breakend circos and return it as inline SVG markup.
#
# Returns list(svg, n_links, n_genes, chroms, interactive, reason). `svg` is NULL when
# there is nothing to draw or too much of it, with `reason` saying which. `interactive`
# is FALSE when the arcs could not all be tagged — the plot is still returned, but the
# caller must not claim it highlights, because it would be highlighting the wrong arcs.
draw_bnd_circos = function(links, genes = NULL, cytobands, chrom_lengths, chromosomes,
                           width = 7.2, height = 7.2) {

  fail = function(reason) list(svg = NULL, n_links = 0L, n_genes = 0L,
                               chroms = character(0), interactive = FALSE,
                               reason = reason)

  if (is.null(links) || nrow(links) == 0)
    return(fail("No breakends with two mapped loci — nothing to draw arcs between."))
  if (nrow(links) > BND_CIRCOS_MAX_LINKS)
    return(fail(sprintf(paste("%s breakend arcs is more than this plot can show",
                              "(limit %s). Select a gene panel above for a filtered",
                              "report, or read the table below."),
                        nrow(links), BND_CIRCOS_MAX_LINKS)))

  b = copy(as.data.table(links))
  chroms = chromosomes[chromosomes %in% unique(c(b$chrom_a, b$chrom_b))]
  chroms = chroms[chroms %in% unique(cytobands$chrom)]
  b = b[chrom_a %in% chroms & chrom_b %in% chroms]
  if (nrow(b) == 0 || length(chroms) == 0)
    return(fail("No breakends on the chromosomes this report plots."))
  chroms = chromosomes[chromosomes %in% unique(c(b$chrom_a, b$chrom_b))]

  cyto_filt = as.data.frame(cytobands[cytobands$chrom %in% chroms, ])
  lens = chrom_lengths[chroms]
  # A position past the end of its sector makes circos.link() error out, which would
  # break the row-to-arc correspondence; clamp instead.
  b[, `:=`(pos_a = pmin(pmax(as.numeric(pos_a), 1), lens[chrom_a]),
           pos_b = pmin(pmax(as.numeric(pos_b), 1), lens[chrom_b]))]

  g = if (is.null(genes)) data.table() else as.data.table(genes)[chrom %in% chroms]
  if (nrow(g) > 0) setorder(g, chrom, start, gene)

  n_link = nrow(b)
  n_gene = nrow(g)
  n_chr = length(chroms)

  link_sent  = .bnd_sentinel("link",  seq_len(n_link))
  chrom_sent = .bnd_sentinel("chrom", seq_len(n_chr))
  gene_sent  = if (n_gene > 0) .bnd_sentinel("gene", seq_len(n_gene)) else character(0)
  line_sent  = if (n_gene > 0) .bnd_sentinel("line", seq_len(n_gene)) else character(0)
  body_sent  = if (n_gene > 0) .bnd_sentinel("body", seq_len(n_gene)) else character(0)

  # A gene span is invisible at this scale — MYC is 7 kb against a 250 Mb sector — so
  # bodies get a floor of about half a degree of arc. The footnote says they are markers,
  # not spans drawn to scale.
  min_w = sum(as.numeric(lens)) / 360

  s = svglite::svgstring(width = width, height = height, bg = "transparent")
  drew = tryCatch({
    circos.clear()
    circos.par(start.degree = 90,
               gap.degree   = c(rep(2, n_chr - 1), 6),
               cell.padding = c(0, 0, 0, 0),
               track.margin = c(0.006, 0.006))
    circos.initializeWithIdeogram(cyto_filt, chromosome.index = chroms, plotType = NULL)

    # Outside labels have to be laid out before the ideogram they point at.
    if (n_gene > 0) {
      circos.genomicLabels(as.data.frame(g[, .(chrom, start, end, gene)]),
                           labels.column     = 4,
                           side              = "outside",
                           col               = gene_sent,
                           line_col          = line_sent,
                           cex               = 0.55,
                           connection_height = mm_h(4))
    }

    circos.genomicIdeogram(cyto_filt)
    ideo_track = get.all.track.index()[length(get.all.track.index())]

    # Gene bodies are drawn onto the ideogram track itself, so a body sits on the band
    # it overlaps rather than in a ring of its own.
    if (n_gene > 0) {
      for (k in seq_len(n_gene)) {
        xlim = get.cell.meta.data("xlim", sector.index = g$chrom[k],
                                  track.index = ideo_track)
        xl = as.numeric(g$start[k]); xr = as.numeric(g$end[k])
        if (xr - xl < min_w) {
          mid = (xl + xr) / 2
          xl = mid - min_w / 2
          xr = mid + min_w / 2
        }
        xl = max(xl, xlim[1]); xr = min(xr, xlim[2])
        if (xr <= xl) next
        circos.rect(xl, 0, xr, 1, col = body_sent[k], border = NA,
                    sector.index = g$chrom[k], track.index = ideo_track)
      }
    }

    # Chromosome names: plotType = NULL suppressed circlize's own, which is what lets
    # each one carry its own sentinel and so be dimmed when its chromosome goes quiet.
    circos.track(ylim = c(0, 1), track.height = 0.07, bg.border = NA,
                 panel.fun = function(x, y) {
                   ch = CELL_META$sector.index
                   circos.text(CELL_META$xcenter, 0.4, sub("^chr", "", ch),
                               facing = "inside", niceFacing = TRUE, cex = 0.62,
                               col = chrom_sent[match(ch, chroms)])
                 })

    for (i in seq_len(n_link)) {
      circos.link(b$chrom_a[i], b$pos_a[i], b$chrom_b[i], b$pos_b[i],
                  col = link_sent[i], lwd = 1)
    }
    TRUE
  }, error = function(e) {
    message("Breakend circos failed: ", conditionMessage(e))
    FALSE
  })

  try(circos.clear(), silent = TRUE)
  dev.off()
  if (!drew) return(fail("The breakend circos could not be drawn."))

  registry = list()
  for (i in seq_len(n_link)) {
    # Both contigs travel with the arc so the client can dim a chromosome label once
    # nothing visible touches it, without a second lookup table.
    registry[[link_sent[i]]] = list(
      family = "link",
      attrs  = sprintf('class="bnd-link" data-svid="%s" data-chrom-a="%s" data-chrom-b="%s"',
                       .attr(b$id[i]), .attr(b$chrom_a[i]), .attr(b$chrom_b[i])))
  }
  for (i in seq_len(n_chr)) {
    registry[[chrom_sent[i]]] = list(
      family = "chrom",
      attrs  = sprintf('class="bnd-chrom-label" data-chrom="%s"', .attr(chroms[i])))
  }
  for (k in seq_len(n_gene)) {
    ga = sprintf('data-gene="%s" data-panels="%s"', .attr(g$gene[k]), .attr(g$panels[k]))
    registry[[gene_sent[k]]] = list(
      family = "gene", expect_text = as.character(g$gene[k]),
      attrs  = paste0('class="bnd-gene-label" ', ga))
    registry[[line_sent[k]]] = list(
      family = "line",
      attrs  = sprintf('class="bnd-gene-line" data-gene="%s"', .attr(g$gene[k])))
    registry[[body_sent[k]]] = list(
      family = "body",
      attrs  = paste0('class="bnd-gene-body" ', ga))
  }

  tagged = .tag_bnd_svg(s(), registry)

  # Every arc must have been tagged exactly once. A plot that highlights the wrong arc is
  # worse than one that does not highlight at all, so a mismatch downgrades the whole
  # thing to a static figure rather than being papered over.
  link_seen = tagged$seen[tagged$seen %in% link_sent]
  interactive = length(link_seen) == n_link && !any(duplicated(link_seen))
  if (!interactive) {
    message("Breakend circos: tagged ", length(unique(link_seen)), " of ", n_link,
            " arcs; highlighting disabled.")
  }
  if (length(tagged$mismatched) > 0) {
    message("Breakend circos: ", length(tagged$mismatched),
            " gene labels did not match their expected symbol and were left untagged.")
  }

  list(svg         = tagged$svg,
       n_links     = n_link,
       n_genes     = sum(tagged$seen %in% gene_sent),
       chroms      = chroms,
       interactive = interactive,
       reason      = NULL)
}
