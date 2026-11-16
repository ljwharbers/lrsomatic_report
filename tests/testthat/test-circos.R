# Genome circos driver labels: the label rows and loci (R/key_findings.R) and the anchors
# draw_circos() reads off the open svglite device (R/circos.R).

test_that("label rows keep tier 1-2 lines and split a joined row into its members", {
  f = data.table(gene = c("EWSR1–ATF1", "CDKN2A", "MYC", "TP53"),
                 genes = c("EWSR1,ATF1", "CDKN2A", "MYC", "TP53"),
                 panels = c("sarcoma", "sarcoma", "lymphoid", ""),
                 tier = c(1L, 2L, 3L, 1L), type = c("SV", "CNA", "SV", "small variant"),
                 alteration = c("translocation", "deletion", "near MYC", "p.R175H"))
  r = circos_driver_rows(f)
  expect_equal(r$label, c("EWSR1", "ATF1", "CDKN2A", "TP53"))
  expect_equal(r$group[1:2], c("EWSR1–ATF1", "EWSR1–ATF1"))
  # Non-panel rows stay, for a typed Custom… list
  expect_equal(r$panels[4], "")
  expect_equal(nrow(circos_driver_rows(NULL)), 0)
  expect_equal(nrow(circos_driver_rows(f[tier == 3L])), 0)
})

test_that("label loci come from the lists, then CN gene rows, then small variants", {
  panels = list(a = load_gene_panel(file.path(repo_root, "tests", "fixtures", "overview",
                                              "panel_a.hg38.tsv"), "hg38"))
  panels = yaml::yaml.load(yaml::as.yaml(panels))
  sd = list(cna = list(genes = data.table(gene = c("CDKN2A", "TP53"), chrom = "chr9",
                                          start = c(21967752L, 1L), end = c(21995301L, 2L))))
  vt = data.table(symbol = c("KMT2D", "TP53"), chrom = c("chr12", "chr1"),
                  pos = c(49020000L, 5L))
  l = circos_label_loci(c("tp53", "CDKN2A", "KMT2D", "NOPOS"), panels, sd, vt)
  setkey(l, gene)
  # The list's interval wins over the CN row and the variant for TP53, matched case-insensitively
  expect_equal(l["tp53"]$chrom, "chr17")
  expect_equal(l["tp53"]$pos, (7661779 + 7687538) / 2)
  expect_equal(l["CDKN2A"]$chrom, "chr9")
  expect_equal(l["KMT2D"]$pos, 49020000)
  expect_false("NOPOS" %in% l$gene)
  # A NOTCH1 row is `snv`-scoped, without SV intervals: the all-scope coordinates still place it
  expect_equal(circos_label_loci("NOTCH1", panels)$chrom, "chr9")
})

test_that("draw_circos anchors sit on their chromosome, outside the names, in SVG units", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("svglite")
  suppressPackageStartupMessages(library(circlize))
  source(file.path(repo_root, "R/circos.R"))
  cyto = load_cytobands("hg38", file.path(repo_root, "assets"))
  lens = load_chrom_lengths("hg38", file.path(repo_root, "assets"))
  chroms = c("chr1", "chr2", "chr3")
  out = tempfile(fileext = ".svg")
  loci = data.table(gene = c("START1", "MID2", "BAD"), chrom = c("chr1", "chr2", "chrZ"),
                    pos = c(1, lens[["chr2"]] / 2, 5))
  res = draw_circos(cytobands = cyto, chrom_lengths = lens, chromosomes = chroms,
                    output_path = out, label_loci = loci)
  expect_true(file.exists(out))
  a = res$anchors
  g = res$geometry
  # An unplotted contig has no anchor
  expect_equal(a$gene, c("START1", "MID2"))
  svg = readLines(out, n = 3, warn = FALSE)
  expect_true(any(grepl(sprintf("viewBox='0 0 %.2f %.2f'", g$width, g$height), svg)))
  # chr1 starts at 12 o'clock (start.degree 90): the anchor is straight above the centre
  expect_lt(abs(a$x[1] - g$cx), 1)
  expect_lt(a$y[1], g$cy)
  r = sqrt((a$x - g$cx)^2 + (a$y - g$cy)^2)
  expect_equal(r, rep(g$r_anchor, 2), tolerance = 0.02)
  expect_gt(g$r_label, g$r_anchor)
  expect_gt(g$r_max, g$r_label)
  expect_true(all(a$x > 0 & a$x < g$width & a$y > 0 & a$y < g$height))
  # The theta convention JS relies on: (cx + r cos, cy - r sin)
  th = a$theta * pi / 180
  expect_equal(a$x, g$cx + g$r_anchor * cos(th), tolerance = 0.02)
  expect_equal(a$y, g$cy - g$r_anchor * sin(th), tolerance = 0.02)

  # Without loci: default canvas, no anchors
  res0 = draw_circos(cytobands = cyto, chrom_lengths = lens, chromosomes = chroms,
                     output_path = tempfile(fileext = ".svg"))
  expect_null(res0$anchors)
})
