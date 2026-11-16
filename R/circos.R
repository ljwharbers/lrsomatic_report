suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
})

# Colour palettes — keep in sync with --circos-* in assets/styles/report.scss

# SBS-6 SNV palette (SigProfiler/COSMIC standard, softened slightly toward report ink/paper).
# No longer drawn on the genome circos (that track is a density now, see Track 1); kept as
# the report's one definition of the SBS-6 colours, mirrored by --circos-snv-* in the SCSS
SNV_COLOURS = c(
  "C>A" = "#2EBAED",
  "C>G" = "#1b1e22",
  "C>T" = "#b3402f",
  "T>A" = "#c7c2b8",
  "T>C" = "#ADCC54",
  "T>G" = "#F0D0CE"
)

# SV_COLOURS and SV_YPOS are defined in R/parse_severus.R, next to the code that stamps them
# onto each row — see the note there on why there is only one copy.

# CNV colours — tied to the report spine
CNV_COLOURS = c(
  major = "#b3402f",  # brick = "more"
  minor = "#0d5c75",  # teal  = "less"
  total = "#1b1e22"   # ink
)

# BND/translocation link colour
BND_COLOUR = "#8a5fa3"

# Mutation-density track fill: warm grey, deliberately outside the SV/CN/SBS hues. The
# legend swatch in templates/sections/_circos.qmd uses the same literal (no --circos-* var)
DENSITY_COLOUR = "#6b655c"

# Canvas half-width when driver labels are drawn: the ring outside the chromosome names that
# assets/js/circos_drivers.js writes gene labels into (1 = circlize's default, no ring)
CIRCOS_LABEL_CANVAS = 1.34

# Draw a circos plot to output_path (SVG or PNG by extension).
#
# label_loci  data.table(gene, chrom, pos) or NULL: genes the browser may label (one row per
#             gene; circos_label_loci() in R/key_findings.R). When non-empty the canvas widens
#             to CIRCOS_LABEL_CANVAS to leave a label ring, and every placeable gene gets an
#             anchor in SVG user units, read off the open device with grconvertX/Y (svglite's
#             device units are its viewBox units, y already pointing down).
# Returns list(path, anchors = data.table(gene, chrom, pos, theta, x, y) | NULL,
#              geometry = list(width, height, cx, cy, r_anchor, r_label, r_max) | NULL)
#   theta  degrees, counter-clockwise from 3 o'clock (circlize's convention; in SVG space a
#          point at radius r is (cx + r cos theta, cy - r sin theta))
#   r_anchor  just outside the chromosome-name track; r_label where labels start; r_max the
#             canvas edge
draw_circos = function(snv_data = NULL,
                       sv_nontrans = NULL,
                       sv_trans = NULL,
                       cnv_data = NULL,
                       cytobands,
                       chrom_lengths,
                       chromosomes,
                       output_path,
                       label_loci = NULL) {

  # Filter cytobands and lengths to displayed chromosomes
  cyto_filt = cytobands[cytobands$chrom %in% chromosomes, ]
  lens_filt  = chrom_lengths[names(chrom_lengths) %in% chromosomes]
  lens_filt  = lens_filt[chromosomes[chromosomes %in% names(lens_filt)]]

  # Prepare SNV density (bin_snv_density() lives in R/parse_smallvariants.R, circlize-free)
  if (!is.null(snv_data) && nrow(snv_data) > 0) {
    snv = as.data.table(snv_data)[nchar(ref) == 1 & nchar(alt) == 1 & chrom %in% chromosomes]
  } else {
    snv = data.table(chrom = character(), pos = integer(), ref = character(), alt = character())
  }
  dens = bin_snv_density(snv, lens_filt, CIRCOS_SNV_BIN)
  dens_cap = circos_density_cap(dens$per_mb)

  # Prepare SV (non-BND) data
  if (!is.null(sv_nontrans) && nrow(sv_nontrans) > 0) {
    sv_nt = as.data.table(sv_nontrans)[chrom %in% chromosomes]
  } else {
    sv_nt = data.table(chrom = character(), pos = integer(), end = integer(),
                       svtype = character(), circos_pos = numeric(), circos_col = character())
  }

  # Prepare translocation (BND) data
  if (!is.null(sv_trans) && nrow(sv_trans) > 0) {
    sv_tr = as.data.table(sv_trans)[chrom %in% chromosomes & chrom2 %in% chromosomes]
  } else {
    sv_tr = data.table(chrom = character(), pos = integer(),
                       chrom2 = character(), pos2 = integer())
  }

  # Prepare CNV data
  if (!is.null(cnv_data) && nrow(cnv_data) > 0) {
    cnv = as.data.table(cnv_data)[chr %in% chromosomes]
    cnv = cnv[order(chr, startpos)]
  } else {
    cnv = data.table(chr = character(), startpos = integer(), endpos = integer(),
                     major_cn = numeric(), minor_cn = numeric(), total_cn = numeric())
  }

  # Driver-label loci on the plotted chromosomes, inside the chromosome
  lab = if (!is.null(label_loci) && nrow(label_loci) > 0)
          as.data.table(label_loci)[chrom %in% names(lens_filt) & !is.na(pos)] else NULL
  if (!is.null(lab) && nrow(lab) > 0) {
    lab = lab[!duplicated(gene)]
    lab[, pos := pmin(pmax(as.numeric(pos), 0), as.numeric(lens_filt[chrom]))]
  } else lab = NULL
  canvas = if (is.null(lab)) 1 else CIRCOS_LABEL_CANVAS

  # Open device
  ext = tolower(tools::file_ext(output_path))
  if (ext == "svg") {
    svglite::svglite(output_path, width = 8, height = 8)
  } else {
    png(output_path, width = 2400, height = 2400, res = 300)
  }

  plot.new()
  circos.clear()

  n_chr = length(chromosomes)
  gap_degrees = c(rep(1.5, n_chr - 1), 7)

  # Single quiet track surface; colours match the report.scss border/surface tokens
  track_bg     = "#fbfaf6"
  track_border = "#e4e0d6"

  circos.par(
    "start.degree" = 90,
    "gap.degree"   = gap_degrees,
    "track.margin" = c(0.006, 0.006),
    "cell.padding" = c(0, 0, 0, 0),
    # Arc resolution (default 500 segments per full circle). Every arc in the SVG is a
    # flattened point list, so this scales the size of every track background, cytoband
    # and CN rectangle; at 200 the chord error on the 680 px plate is ~0.04 px
    "unit.circle.segments" = 200,
    "canvas.xlim" = c(-canvas, canvas),
    "canvas.ylim" = c(-canvas, canvas)
  )

  # Build cytobands list as expected by circos.initializeWithIdeogram
  cyto_list = list(
    df          = cyto_filt,
    chromosome  = chromosomes[chromosomes %in% unique(cyto_filt$chrom)],
    chr.len     = lens_filt
  )

  circos.initializeWithIdeogram(cyto_list$df,
                                chromosome.index = cyto_list$chromosome,
                                plotType         = c("ideogram", "labels"),
                                labels.cex       = 0.8)

  # ---- Track 1: SNV density (SNVs per CIRCOS_SNV_BIN window) ---------------
  # One filled area per chromosome, a vertex per window, so the SVG size no longer scales
  # with mutation load (one point per SNV was 5-9 MB). Clipped at `dens_cap`: a single
  # kataegis window would otherwise flatten every other bin to the baseline.
  circos.trackPlotRegion(
    factors      = chromosomes,
    ylim         = c(0, dens_cap),
    bg.border    = track_border,
    bg.col       = track_bg,
    track.height = 0.16,
    panel.fun    = function(region, value, ...) {
      chr = get.cell.meta.data("sector.index")
      d = dens[chrom == chr]
      if (nrow(d) == 0 || all(d$per_mb == 0)) return(invisible(NULL))
      # One vertex per window, at its midpoint, held flat to both chromosome ends. Interior
      # vertices of a run of equal (capped) heights are dropped: they add nothing to the
      # outline, and in a low-burden sample most windows are a run of zeros
      y = pmin(d$per_mb, dens_cap)
      x = c(0, (d$start + d$end) / 2, d$end[nrow(d)])
      y = c(y[1], y, y[length(y)])
      k = length(y)
      keep = c(TRUE, y[-1] != y[-k]) | c(y[-k] != y[-1], TRUE)
      x = x[keep]; y = y[keep]
      # Closed along the baseline; every 10th vertex is enough for it to follow the arc
      # (a 10 Mb chord is <0.05 px off the curve on the 680 px plate)
      bx = rev(unique(c(x[seq(1, length(x), by = 10)], x[length(x)])))
      pol = circlize(c(x, bx), c(y, rep(0, length(bx))))
      # A plain polygon in canvas coordinates: circos.lines(area = TRUE) subdivides every
      # segment to follow the curvature, which tripled this track's SVG weight
      th = pol[, "theta"] * pi / 180
      polygon(pol[, "rou"] * cos(th), pol[, "rou"] * sin(th),
              col = adjustcolor(DENSITY_COLOUR, alpha.f = 0.55),
              border = DENSITY_COLOUR, lwd = 0.4)
    }
  )

  dens_track_index = get.current.track.index()
  if (nrow(snv) > 0) tryCatch(
    circos.yaxis(
      side              = "left",
      at                = c(0, dens_cap),
      labels            = c("0", paste0(formatC(dens_cap, format = "fg", digits = 3), "/Mb")),
      track.index       = dens_track_index,
      sector.index      = chromosomes[1],
      labels.niceFacing = TRUE,
      labels.cex        = 0.40
    ),
    error = function(e) NULL
  )

  # ---- Track 2: Non-BND SVs (DEL/DUP/INV/INS as horizontal segments) -----
  circos.trackPlotRegion(
    factors      = chromosomes,
    ylim         = c(0, 1),
    bg.border    = track_border,
    bg.col       = track_bg,
    track.height = 0.10,
    panel.fun    = function(region, value, ...) {
      chr = get.cell.meta.data("sector.index")
      sub_sv = sv_nt[chrom == chr & !is.na(circos_pos)]
      if (nrow(sub_sv) == 0) return(invisible(NULL))
      for (i in seq_len(nrow(sub_sv))) {
        x1 = sub_sv$pos[i]
        x2 = if (!is.na(sub_sv$end[i]) && sub_sv$end[i] > x1) sub_sv$end[i] else x1 + 1L
        circos.segments(
          x0  = x1, x1  = x2,
          y0  = sub_sv$circos_pos[i], y1 = sub_sv$circos_pos[i],
          col = sub_sv$circos_col[i],
          lwd = 2.5
        )
      }
    }
  )

  # Ask circlize which track that was rather than counting.
  # circos.initializeWithIdeogram() creates a track per plotType group *before* any of
  # ours -- one for axis/labels and one for the ideogram -- so the SV track is index 4,
  # not 3. A wrong track.index does not error, it silently draws on the neighbouring
  # track, which is how these labels ended up against the SNV dots.
  sv_track_index = get.current.track.index()

  # Y-axis labels for SV track. `at`/`labels` come from SV_YPOS (R/parse_severus.R) — the
  # same object severus_circos_tracks() positioned the segments from, not a second copy of
  # its values, so the two cannot drift apart.
  sv_at = sort(SV_YPOS)
  tryCatch(
    circos.yaxis(
      side              = "left",
      at                = unname(sv_at),
      labels            = names(sv_at),
      track.index       = sv_track_index,
      sector.index      = chromosomes[1],
      labels.niceFacing = TRUE,
      labels.cex        = 0.45
    ),
    error = function(e) NULL
  )

  # ---- Track 3: ASCAT copy-number -----------------------------------------
  circos.trackPlotRegion(
    factors      = chromosomes,
    ylim         = c(0, 4),
    bg.border    = track_border,
    bg.col       = track_bg,
    track.height = 0.18,
    panel.fun    = function(region, value, ...) {
      chr = get.cell.meta.data("sector.index")
      sub_cnv = cnv[chr == get.cell.meta.data("sector.index")]
      if (nrow(sub_cnv) == 0) return(invisible(NULL))

      xmax = lens_filt[chr]
      if (!is.na(xmax)) {
        for (y_ref in c(1, 2, 3)) {
          circos.lines(c(0, xmax), c(y_ref, y_ref),
                       col = track_border, lwd = 0.4, lty = "dotted")
        }
      }

      for (i in seq_len(nrow(sub_cnv))) {
        xl = sub_cnv$startpos[i]; xr = sub_cnv$endpos[i]
        maj = sub_cnv$major_cn[i]
        circos.rect(xl, maj + 0.02, xr, maj + 0.12,
                    col = CNV_COLOURS["major"], border = CNV_COLOURS["major"], lwd = 0.05)
        min_cn = sub_cnv$minor_cn[i]
        circos.rect(xl, min_cn - 0.12, xr, min_cn - 0.02,
                    col = CNV_COLOURS["minor"], border = CNV_COLOURS["minor"], lwd = 0.05)
        tot = sub_cnv$total_cn[i]
        circos.rect(xl, tot - 0.03, xr, tot + 0.03,
                    col = CNV_COLOURS["total"], border = CNV_COLOURS["total"], lwd = 0.05)
      }
    }
  )

  # Y-axis labels for the copy-number track. Drawn once, outside panel.fun: in there it
  # ran per sector with a fixed sector.index, so it was re-drawn once for every
  # chromosome carrying ASCAT segments (~23 overlapping copies, which is why it rendered
  # heavier than the SV axis). The empty case keeps its old behaviour deliberately: the
  # track is created unconditionally, but parse_ascat_segments() returns NULL when there is
  # no ASCAT output at all, and a 0..4+ scale against a blank ring reads as "copy number is
  # zero everywhere" rather than "no copy-number data".
  cnv_track_index = get.current.track.index()
  if (nrow(cnv) > 0) tryCatch(
    circos.yaxis(
      side              = "left",
      at                = c(0, 1, 2, 3, 4),
      labels            = c("0", "1", "2", "3", "4+"),
      track.index       = cnv_track_index,
      sector.index      = chromosomes[1],
      labels.niceFacing = TRUE,
      labels.cex        = 0.40
    ),
    error = function(e) NULL
  )

  # ---- Translocation links (BND): one arc per mate-collapsed rearrangement ----
  if (nrow(sv_tr) > 0) {
    for (i in seq_len(nrow(sv_tr))) {
      tryCatch(
        circos.link(
          sector.index1 = sv_tr$chrom[i],  point1 = sv_tr$pos[i],
          sector.index2 = sv_tr$chrom2[i], point2 = sv_tr$pos2[i],
          col = adjustcolor(BND_COLOUR, alpha.f = 0.45),
          lwd = 0.9
        ),
        error = function(e) NULL
      )
    }
  }

  # Driver-label anchors, in device units = SVG user units (svglite). Track 1 is the
  # chromosome-name track circos.initializeWithIdeogram() made, the outermost one
  anchors = NULL; geometry = NULL
  if (!is.null(lab) && ext == "svg") {
    s1 = get.all.sector.index()[1]
    r_names = get.cell.meta.data("cell.top.radius", sector.index = s1, track.index = 1)
    # circlize() takes one sector per call
    pol = do.call(rbind, lapply(seq_len(nrow(lab)), function(i)
      circlize(lab$pos[i], 0, sector.index = lab$chrom[i], track.index = 1)))
    theta = pol[, "theta"]
    dev_x = function(v) grconvertX(v, "user", "device")
    dev_y = function(v) grconvertY(v, "user", "device")
    cx = dev_x(0); cy = dev_y(0)
    unit = dev_x(1) - cx
    r_anchor = (r_names + 0.015) * unit
    th = theta * pi / 180
    anchors = data.table(gene = lab$gene, chrom = lab$chrom, pos = lab$pos,
                         theta = round(theta, 3),
                         x = round(dev_x(r_anchor / unit * cos(th)), 2),
                         y = round(dev_y(r_anchor / unit * sin(th)), 2))
    sz = dev.size("in") * 72
    geometry = list(width = round(sz[1], 2), height = round(sz[2], 2),
                    cx = round(cx, 2), cy = round(cy, 2),
                    r_anchor = round(r_anchor, 2),
                    r_label  = round((r_names + 0.07) * unit, 2),
                    r_max    = round(canvas * unit, 2))
  }

  circos.clear()
  dev.off()
  invisible(list(path = output_path, anchors = anchors, geometry = geometry))
}
