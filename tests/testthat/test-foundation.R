# v2 shared scaffolding: resolve_columns(), count_small_variants(), embed_svg_plot(),
# sv_panel_hit_table(), the section ordering / ctx contract, and the small-variant sidecar.

# ---- resolve_columns ---------------------------------------------------------------

test_that("resolve_columns prefers an exact name, then a normalised one", {
  dt = data.table(Chrom = "chr1", POS = 1L, `sv-type` = "DEL")
  got = resolve_columns(dt, list(chrom = c("chrom", "chr"), pos = "pos",
                                 svtype = c("SVTYPE", "svtype"), vaf = "VAF"))
  expect_equal(got, c(chrom = "Chrom", pos = "POS", svtype = "sv-type", vaf = NA_character_))
})

test_that("resolve_columns maps pos_1 and Pos_1 to their own canonicals", {
  dt = data.table(pos_1 = 1L, Pos_1 = 2L)
  got = resolve_columns(dt, c(a = "Pos_1", b = "pos_1"))
  expect_equal(got, c(a = "Pos_1", b = "pos_1"))
  # Both canonicals reaching the pair only by normalised name still get distinct columns
  got = resolve_columns(dt, list(a = "POS_1", b = "POS_1"))
  expect_equal(unname(sort(got)), sort(c("Pos_1", "pos_1")))
  # An exact claim is never handed to a later canonical's normalised match
  got = resolve_columns(dt, list(first = "posone", x = "POS_1", y = "pos_1"))
  expect_equal(got, c(first = NA_character_, x = "Pos_1", y = "pos_1"))
})

test_that("resolve_columns takes candidates in order and accepts a names vector", {
  got = resolve_columns(c("chr", "chromosome"), list(chrom = c("chromosome", "chr")))
  expect_equal(got, c(chrom = "chromosome"))
  expect_error(resolve_columns(data.table(a = 1), list("a")), "named")
})

# ---- count_small_variants ----------------------------------------------------------

test_that("count_small_variants splits SNVs from indels over unique variants", {
  v = data.table(chrom = c("chr1", "chr1", "chr1", "chr2", "chr3", "chr4"),
                 pos   = c(10L, 10L, 20L, 30L, 40L, 50L),
                 ref   = c("A", "A", "AT", "-", "c", "G"),
                 alt   = c("G", "G", "A", "T", "t", "N"),
                 symbol = c("X", "Y", "X", "X", "X", "X"))
  got = count_small_variants(v)
  # chr1:10 A>G is one variant over two genes; "-"/"T" (VEP dash-form insertion) is an
  # indel; lower-case bases count; N is not an ACGT base
  expect_equal(got, list(total = 5L, snv = 2L, indel = 3L))
})

test_that("count_small_variants is NA without a VEP file and zero for an empty one", {
  expect_equal(count_small_variants(NULL),
               list(total = NA_integer_, snv = NA_integer_, indel = NA_integer_))
  empty = data.table(chrom = character(), pos = integer(), ref = character(), alt = character())
  expect_equal(count_small_variants(empty), list(total = 0L, snv = 0L, indel = 0L))
})

# ---- embed_svg_plot ----------------------------------------------------------------

test_that("embed_svg_plot inlines a ggplot as a base64 SVG <img>", {
  skip_if_not_installed("svglite")
  p = ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  tag = embed_svg_plot(p, width = 4, height = 3, alt = "three points")
  html = as.character(tag)
  expect_match(html, "^<img ")
  expect_match(html, 'src="data:image/svg\\+xml;base64,')
  expect_match(html, 'alt="three points"')
  expect_match(html, 'class="report-figure"')
  b64 = sub('.*base64,([^"]+)".*', "\\1", html)
  svg = rawToChar(base64enc::base64decode(b64))
  expect_match(svg, "<svg")
  expect_null(embed_svg_plot(NULL, 4, 3))
})

test_that("embed_svg_plot accepts a base-graphics function and closes its device", {
  skip_if_not_installed("svglite")
  n_dev = length(grDevices::dev.list())
  tag = embed_svg_plot(function() plot(1:3), width = 3, height = 3, class = "x")
  expect_match(as.character(tag), 'class="x"')
  expect_equal(length(grDevices::dev.list()), n_dev)
  # A plot that errors still closes the device
  expect_error(embed_svg_plot(function() stop("boom"), 3, 3), "boom")
  expect_equal(length(grDevices::dev.list()), n_dev)
})

# ---- sv_panel_hit_table ------------------------------------------------------------

hit_sv = function() data.table(
  id      = c("bnd1", "del1", "bnd2"),
  svtype  = c("BND", "DEL", "BND"),
  chrom_a = c("chr8", "chr10", "chr3"),
  pos_a   = c(1000000L, 5000000L, 100L),
  chrom_b = c("chr14", "chr10", "chr3"),
  pos_b   = c(2000000L, 5100000L, 900L),
  gene_a  = c("MYC,PVT1", NA, "BCL6"),
  gene_b  = c("IGH", NA, NA)
)
hit_coord_panel = function(genes, chrom, start, end, name = "p") list(
  name = name, has_coords = TRUE, genes = unique(genes), genes_sv = unique(genes),
  chrom = chrom, start = as.integer(start), end = as.integer(end), interval_gene = genes)

test_that("sv_panel_hit_table is the long form of sv_panel_hits' labels", {
  sv = hit_sv()
  coord = hit_coord_panel(c("INSIDE", "NEAR", "BSIDE", "SPANNED", "TWICE", "TWICE"),
                          c("chr8", "chr8", "chr14", "chr10", "chr3", "chr3"),
                          c(999000, 1300000, 2000500, 5040000, 50, 5000),
                          c(1001000, 1310000, 2001000, 5060000, 150, 6000), name = "coord")
  sym = list(name = "sym", has_coords = FALSE, genes = c("PVT1", "MYC", "IGH", "BCL6"))

  ht = sv_panel_hit_table(sv, coord)
  expect_equal(names(ht), c("row", "id", "gene", "side", "gap", "direct", "panel", "match"))
  expect_equal(ht[gene == "INSIDE", .(row, id, side, gap, direct)],
               data.table(row = 1L, id = "bnd1", side = "A", gap = 0, direct = TRUE))
  expect_equal(ht[gene == "NEAR"]$gap, 300000)
  expect_equal(ht[gene == "BSIDE"]$side, "B")
  expect_equal(ht[gene == "SPANNED"]$side, "span")
  # One gene with two loci near one breakend: two rows, two gaps
  expect_equal(sort(ht[gene == "TWICE" & side == "A"]$gap), c(0, 4900))
  expect_true(all(ht$panel == "coord"))

  hs = sv_panel_hit_table(sv, sym)
  expect_true(all(hs$match == "symbol" & hs$direct & hs$gap == 0))
  # Symbol-mode genes are upper-cased VEP symbols in cell order
  expect_equal(hs[row == 1 & side == "A"]$gene, c("MYC", "PVT1"))

  # Labels rebuilt from the long form, for one panel, several, and mixed modes
  for (ps in list(coord, list(coord = coord), sym, list(coord = coord, sym = sym),
                  list(sym = sym, coord = coord))) {
    labs = sv_panel_hits(sv, ps)
    ht   = sv_panel_hit_table(sv, ps)
    tagged = length(.as_panel_list(ps)) > 1
    rebuilt = rep("", nrow(sv))
    lab = .sv_hit_labels(ht, tagged)
    for (r in unique(lab$row)) rebuilt[r] = paste(unique(lab[row == r]$label), collapse = ", ")
    expect_identical(rebuilt, labs)
  }
  expect_equal(sv_panel_hits(sv, sym)[1], "MYC,PVT1 (A, direct), IGH (B, direct)")
  expect_equal(sv_panel_hits(sv, list(coord = coord, sym = sym))[2],
               "SPANNED (span, direct) [coord]")
})

test_that("sv_panel_hit_table is empty-safe", {
  ht = sv_panel_hit_table(hit_sv(), NULL)
  expect_equal(nrow(ht), 0)
  expect_equal(names(ht), c("row", "id", "gene", "side", "gap", "direct", "panel", "match"))
  expect_equal(nrow(sv_panel_hit_table(data.table(), hit_coord_panel("X", "chr1", 1, 2))), 0)
  expect_equal(nrow(sv_panel_hit_table(hit_sv(), hit_coord_panel("X", "chr1", 1, 2))), 0)
})

# ---- section ordering and ctx ------------------------------------------------------

test_that("section modules are sourced and registered in SECTION_ORDER", {
  files = section_files(file.path(repo_root, "R/sections"))
  ids = sub("\\.R$", "", basename(files))
  expect_equal(ids[seq_along(SECTION_ORDER)][ids[seq_along(SECTION_ORDER)] %in% SECTION_ORDER],
               SECTION_ORDER[SECTION_ORDER %in% ids])
  expect_equal(names(SECTIONS)[seq_along(SECTION_ORDER)], SECTION_ORDER)
})

test_that("an unlisted section file is sourced after the ordered ones", {
  d = tempfile("sections")
  dir.create(d)
  for (f in c("aaa.R", "whatshap.R", "sv.R", "zzz.R")) file.create(file.path(d, f))
  expect_equal(basename(section_files(d)), c("sv.R", "whatshap.R", "aaa.R", "zzz.R"))
})

test_that("ctx is passed only to a parse() that declares it", {
  two   = list(parse = function(inputs, section_data) list(n = length(section_data)))
  three = list(parse = function(inputs, section_data, ctx) ctx$reference)
  expect_equal(run_section_parse(two, list(), list(a = 1), list(reference = "hg38")),
               list(n = 1L))
  expect_equal(run_section_parse(three, list(), list(), list(reference = "hg38")), "hg38")
})

test_that("every stub section parses to NULL harmlessly", {
  for (id in c("cna", "shatterseek", "ecdna", "reconplot", "signatures")) {
    s = SECTIONS[[id]]
    expect_false(is.null(s), info = id)
    expect_null(run_section_parse(s, s$locate(tempdir(), "S"), list(), list()), info = id)
  }
})

test_that("build_overview's stub returns the findings contract", {
  o = build_overview(NULL, list(), list(), list(purity = NA, ploidy = NA),
                     list(tmb = NA), list(), list())
  expect_equal(names(o), c("findings", "biomarkers", "per_panel_n"))
  expect_equal(names(o$findings), c("gene", "type", "alteration", "evidence", "tier",
                                    "role", "panels", "vaf", "source_id"))
  expect_type(o$findings$tier, "integer")
})

# ---- small-variant embedding and sidecar --------------------------------------------

test_that("slim_variant_table keeps attributes and records the full row count", {
  vt = data.table(symbol = c("A", "B"), impact = c("HIGH", "LOW"))
  setattr(vt, "vep_fields", c(sift = "SIFT"))
  out = slim_variant_table(vt, list(), "slim")
  expect_equal(attr(out, "vep_fields"), c(sift = "SIFT"))
  expect_equal(attr(out, "n_full"), 2L)
  # The caller's full table is not tagged by reference
  expect_null(attr(vt, "n_full"))
  expect_null(slim_variant_table(NULL, list()))
})

test_that("write_full_variant_table writes a gzipped TSV and drops list columns", {
  vt = data.table(symbol = c("A", "B"), pos = 1:2)
  vt[, extra := list(list(1, 2))]
  f = tempfile(fileext = ".tsv.gz")
  write_full_variant_table(vt, f)
  back = fread(f)
  expect_equal(names(back), c("symbol", "pos"))
  expect_equal(back$pos, 1:2)
  expect_null(write_full_variant_table(NULL, tempfile()))
})
