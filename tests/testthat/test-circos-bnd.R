assets_dir = file.path(repo_root, "assets")

cyto_hg38 = function() load_cytobands("hg38", assets_dir)
lens_hg38 = function() load_chrom_lengths("hg38", assets_dir)
chroms_hg38 = function() {
  ch = chromosomes_for_sex("female")
  ch[ch %in% unique(cyto_hg38()$chrom)]
}

# Three drawable rearrangements plus two rows that must never reach the plot: a DEL
# (not a breakend) and a single breakend (no partner locus).
sv_fixture = function() {
  data.table(
    id      = c("severus_BND0_1", "severus_BND1_1", "severus_BND2_1",
                "severus_5", "severus_BND9_1"),
    svclass = c("translocation", "translocation", "intra-chr breakend",
                "DEL", "single breakend"),
    chrom_a = c("chr8", "chr14", "chr18", "chr1", "chr3"),
    pos_a   = c(127736000L, 105586000L, 63123500L, 1000000L, 500000L),
    chrom_b = c("chr14", "chr18", "chr18", "chr1", NA_character_),
    pos_b   = c(105590000L, 63200000L, 70000000L, 1050000L, NA_integer_))
}

n_class = function(svg, cls) {
  m = gregexpr(paste0('class="', cls, '"'), svg, fixed = TRUE)[[1]]
  if (identical(as.integer(m[1]), -1L)) 0L else length(m)
}

# ---- bnd_links -----------------------------------------------------------

test_that("bnd_links keeps only breakends with two mapped loci", {
  bl = bnd_links(sv_fixture(), chroms_hg38())
  expect_equal(nrow(bl), 3)
  expect_setequal(bl$id, c("severus_BND0_1", "severus_BND1_1", "severus_BND2_1"))
  expect_false("severus_5" %in% bl$id)        # a DEL is not a breakend
  expect_false("severus_BND9_1" %in% bl$id)   # a single breakend has no partner
})

test_that("bnd_links drops breakends off the plotted chromosome set", {
  bl = bnd_links(sv_fixture(), c("chr8", "chr14"))
  expect_equal(bl$id, "severus_BND0_1")
})

test_that("bnd_links tolerates an empty or malformed table", {
  expect_equal(nrow(bnd_links(NULL, chroms_hg38())), 0)
  expect_equal(nrow(bnd_links(data.table(), chroms_hg38())), 0)
  expect_equal(nrow(bnd_links(data.table(id = "x"), chroms_hg38())), 0)
})

# ---- bnd_panel_genes -----------------------------------------------------

test_that("bnd_panel_genes matches on the same window as the panel_hit column", {
  panels = load_all_gene_panels(assets_dir, "hg38")
  skip_if(is.null(panels$lymphoid) || !isTRUE(panels$lymphoid$has_coords))

  bl = bnd_links(sv_fixture(), chroms_hg38())
  g  = bnd_panel_genes(bl, panels)

  # MYC (chr8:127,735,434) is 1.5 kb from the chr8 breakend; BCL2 (chr18:63,123,346)
  # sits under the chr18 one.
  expect_true(all(c("MYC", "BCL2") %in% g$gene))
  expect_true(all(g$panels == "lymphoid"))
  expect_equal(anyDuplicated(g[, .(chrom, start, end, gene)]), 0L)

  # Half the window away still hits; well outside it does not.
  near = data.table(id = "a", svclass = "translocation",
                    chrom_a = "chr8", pos_a = 127735434L - 5e5,
                    chrom_b = "chr14", pos_b = 20000000L)
  far  = copy(near)[, pos_a := 127735434L - 5e6]
  expect_true("MYC" %in% bnd_panel_genes(near, panels)$gene)
  expect_false("MYC" %in% bnd_panel_genes(far, panels)$gene)
})

test_that("bnd_panel_genes returns nothing for a symbol-only panel", {
  symbol_only = list(custom = list(has_coords = FALSE, genes = c("MYC", "BCL2")))
  bl = bnd_links(sv_fixture(), chroms_hg38())
  expect_equal(nrow(bnd_panel_genes(bl, symbol_only)), 0)
  expect_equal(nrow(bnd_panel_genes(bl, list())), 0)
  expect_equal(nrow(bnd_panel_genes(bl, NULL)), 0)
})

# ---- draw_bnd_circos -----------------------------------------------------

test_that("draw_bnd_circos tags every arc and gene, and leaks no sentinel", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("svglite")

  panels = load_all_gene_panels(assets_dir, "hg38")
  bl  = bnd_links(sv_fixture(), chroms_hg38())
  g   = bnd_panel_genes(bl, panels)
  res = draw_bnd_circos(bl, g, cyto_hg38(), lens_hg38(), chroms_hg38())

  expect_true(res$interactive)
  expect_equal(res$n_links, 3L)
  expect_equal(n_class(res$svg, "bnd-link"), 3)

  # Every row's id reached the SVG, and nothing else did.
  svids = unlist(regmatches(res$svg, gregexpr('data-svid="[^"]*"', res$svg)))
  expect_setequal(svids, sprintf('data-svid="%s"', bl$id))

  # Sectors are the touched chromosomes only, in the report's canonical order.
  expect_equal(res$chroms, c("chr8", "chr14", "chr18"))
  expect_equal(n_class(res$svg, "bnd-chrom-label"), 3)

  expect_equal(res$n_genes, nrow(g))
  expect_equal(n_class(res$svg, "bnd-gene-body"), nrow(g))
  expect_equal(n_class(res$svg, "bnd-gene-label"), nrow(g))
  expect_true(n_class(res$svg, "bnd-gene-line") >= nrow(g))
  genes_out = unique(unlist(regmatches(res$svg, gregexpr('data-gene="[^"]*"', res$svg))))
  expect_setequal(genes_out, sprintf('data-gene="%s"', sort(unique(g$gene))))

  # The whole point of the rewrite: no drawing colour survives into the output, and
  # svglite's own global stylesheet no longer targets bare element names.
  expect_false(grepl(BND_SENTINEL_RE, res$svg))
  expect_false(grepl(".svglite ", res$svg, fixed = TRUE))
  expect_true(startsWith(res$svg, "<svg"))
})

test_that("draw_bnd_circos renders without genes", {
  skip_if_not_installed("circlize")
  bl  = bnd_links(sv_fixture(), chroms_hg38())
  res = draw_bnd_circos(bl, NULL, cyto_hg38(), lens_hg38(), chroms_hg38())

  expect_true(res$interactive)
  expect_equal(res$n_genes, 0L)
  expect_equal(n_class(res$svg, "bnd-gene-body"), 0)
  expect_equal(n_class(res$svg, "bnd-gene-label"), 0)
})

test_that("draw_bnd_circos declines rather than drawing nothing or everything", {
  skip_if_not_installed("circlize")
  chroms = chroms_hg38()

  none = draw_bnd_circos(bnd_links(sv_fixture()[0], chroms), NULL,
                         cyto_hg38(), lens_hg38(), chroms)
  expect_null(none$svg)
  expect_match(none$reason, "nothing to draw")

  bl  = bnd_links(sv_fixture(), chroms)
  big = bl[rep(1L, BND_CIRCOS_MAX_LINKS + 1L)]
  over = draw_bnd_circos(big, NULL, cyto_hg38(), lens_hg38(), chroms)
  expect_null(over$svg)
  expect_match(over$reason, as.character(BND_CIRCOS_MAX_LINKS))
})

# ---- .tag_bnd_svg --------------------------------------------------------

test_that(".tag_bnd_svg leaves colours it does not own alone", {
  # #101010 matches the sentinel pattern but is a plausible cytoband grey, so it must
  # survive untouched; only a colour in the registry is rewritten.
  txt = paste(
    "<svg>",
    "<polygon points='0,0' style='fill: #101010;' />",
    "<polyline points='0,0' style='stroke: #100001;' />",
    "</svg>", sep = "\n")
  out = .tag_bnd_svg(txt, list(`#100001` = list(family = "link", attrs = 'class="bnd-link"')))

  expect_true(grepl("#101010", out$svg, fixed = TRUE))
  expect_false(grepl("#100001", out$svg, fixed = TRUE))
  expect_equal(out$seen, "#100001")
  expect_true(grepl('<polyline class="bnd-link"', out$svg, fixed = TRUE))
})

test_that(".tag_bnd_svg refuses a label that is not the one it expected", {
  txt = "<svg>\n<text style='fill: #300001;'>NOTMYC</text>\n</svg>"
  out = .tag_bnd_svg(txt, list(
    `#300001` = list(family = "gene", expect_text = "MYC",
                     attrs = 'class="bnd-gene-label" data-gene="MYC"')))

  expect_equal(out$mismatched, "#300001")
  expect_length(out$seen, 0)
  expect_false(grepl("bnd-gene-label", out$svg, fixed = TRUE))
})
