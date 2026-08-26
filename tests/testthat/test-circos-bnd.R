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
  # The gene track is the superset over every coordinate-carrying panel, so `panels`
  # names all of them — don't pin it to the builtin set, which grows. A gene must
  # come from at least one panel that actually lists it.
  coord_panels = names(panels)[vapply(panels, function(p) isTRUE(p$has_coords), logical(1))]
  for (i in seq_len(nrow(g))) {
    named = strsplit(g$panels[i], ",", fixed = TRUE)[[1]]
    expect_true(length(named) > 0)
    expect_true(all(named %in% coord_panels))
    expect_true(all(vapply(named, function(nm) g$gene[i] %in% panels[[nm]]$genes,
                           logical(1))))
  }
  expect_true("lymphoid" %in% strsplit(g$panels[g$gene == "BCL2"], ",")[[1]])
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

# ---- bnd_circos_data -----------------------------------------------------
#
# The plot itself is drawn by assets/js/bnd_circos.js; what R is responsible for is
# selecting the right records and serialising them correctly, which is what these cover.

test_that("bnd_circos_data emits every drawable arc, over only the chromosomes touched", {
  chroms = chroms_hg38()
  bl  = bnd_links(sv_fixture(), chroms)
  res = bnd_circos_data(bl, bnd_panel_genes(bl, load_all_gene_panels(assets_dir, "hg38")),
                        cyto_hg38(), lens_hg38(), chroms)

  expect_null(res$reason)
  expect_equal(res$n_links, 3L)
  # Only the chromosomes a drawable breakend touches, in reference order.
  expect_equal(res$chroms, c("chr8", "chr14", "chr18"))
  expect_equal(res$data$chromosomes, '["chr8","chr14","chr18"]')
  expect_equal(length(gregexpr('["severus', res$data$links, fixed = TRUE)[[1]]), 3L)

  # Cytobands are restricted to the plotted sectors: the payload is inlined into a
  # self-contained HTML file, so shipping all 1549 hg38 bands would be waste.
  expect_false(grepl('["chr1",', res$data$cytobands, fixed = TRUE))
  expect_true(grepl('["chr8",', res$data$cytobands, fixed = TRUE))
})

test_that("bnd_circos_data emits an empty gene array when no panel carries coordinates", {
  chroms = chroms_hg38()
  res = bnd_circos_data(bnd_links(sv_fixture(), chroms), NULL,
                        cyto_hg38(), lens_hg38(), chroms)
  expect_equal(res$data$genes, "[]")
  expect_equal(res$n_genes, 0L)
})

test_that("bnd_circos_data declines rather than emitting nothing or everything", {
  chroms = chroms_hg38()

  none = bnd_circos_data(bnd_links(sv_fixture()[0], chroms), NULL,
                         cyto_hg38(), lens_hg38(), chroms)
  expect_null(none$data)
  expect_match(none$reason, "nothing to draw")

  big = bnd_links(sv_fixture(), chroms)
  big = rbindlist(rep(list(big), ceiling((BND_CIRCOS_MAX_LINKS + 1) / nrow(big))))
  over = bnd_circos_data(big, NULL, cyto_hg38(), lens_hg38(), chroms)
  expect_null(over$data)
  expect_match(over$reason, "more than this plot can show")
})

test_that("bnd_circos_data clamps a locus past the end of its sector", {
  chroms = chroms_hg38()
  sv = sv_fixture()[1]
  sv$pos_b = 999000000L   # past the end of chr14
  res = bnd_circos_data(bnd_links(sv, chroms), NULL, cyto_hg38(), lens_hg38(), chroms)
  # Clamped rather than dropped: dropping it would break the arc's cross-link to its row.
  expect_equal(res$n_links, 1L)
  expect_true(grepl(format(lens_hg38()[["chr14"]], scientific = FALSE),
                    res$data$links, fixed = TRUE))
})

test_that("bnd_circos_script produces one assignable literal, and nothing when refused", {
  chroms = chroms_hg38()
  res = bnd_circos_data(bnd_links(sv_fixture(), chroms), NULL,
                        cyto_hg38(), lens_hg38(), chroms)
  js = bnd_circos_script(res)
  expect_match(js, "^<script>window\\.BND_DATA=\\{")
  expect_match(js, "\\};</script>$")

  expect_equal(bnd_circos_script(list(data = NULL)), "")
})

# ---- js_vec / js_rows ----------------------------------------------------

test_that("js_vec and js_rows emit valid literals, with NA as null", {
  expect_equal(js_vec(character(0)), "[]")
  expect_equal(js_vec(c("a", NA)), '["a",null]')
  expect_equal(js_vec(c(1, 2.5)), "[1,2.5]")
  # A coordinate must not come out in scientific notation: 2.48956422e+08 is valid JS but
  # reads back as a different number than the base pair it names.
  expect_equal(js_vec(248956422), "[248956422]")

  dt = data.table(a = c("x", "y"), b = c(1L, NA_integer_))
  expect_equal(js_rows(dt, c("a", "b")), '[["x",1],["y",null]]')
  expect_equal(js_rows(data.table(), "a"), "[]")
})

