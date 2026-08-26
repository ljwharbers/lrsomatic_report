test_that("setup loads without error", {
  expect_true(is.function(load_gene_panel))
})

test_that("load_all_gene_panels resolves per-reference builtins to one entry", {
  # Use the real assets dir
  # When testthat runs, getwd() is in tests/testthat, so we need to go up two levels to repo root
  repo_root = dirname(dirname(getwd()))
  for (ref in c("hg38", "t2t")) {
    panels = load_all_gene_panels(file.path(repo_root, "assets"), ref)
    expect_type(panels, "list")
    expect_true("lymphoid" %in% names(panels))
    p = panels[["lymphoid"]]
    expect_true(p$has_coords)
    expect_equal(p$reference, ref)
    expect_true("MYC" %in% p$genes)
    expect_equal(length(p$chrom), length(p$interval_gene))
  }
})

test_that("every builtin panel ships for both references with matching gene sets", {
  # A builtin missing from one reference would be silently unavailable there, and a
  # gene placed on a different chromosome between the two files means a bad liftover.
  repo_root = dirname(dirname(getwd()))
  assets    = file.path(repo_root, "assets")
  hg38      = load_all_gene_panels(assets, "hg38")
  t2t       = load_all_gene_panels(assets, "t2t")
  expect_setequal(names(hg38), names(t2t))
  expect_true(all(c("lymphoid", "sarcoma") %in% names(hg38)))

  for (nm in names(hg38)) {
    a = hg38[[nm]]; b = t2t[[nm]]
    expect_setequal(a$genes, b$genes)
    expect_equal(a$has_coords, b$has_coords)
    if (!isTRUE(a$has_coords)) next
    chrom_a = setNames(a$chrom, a$interval_gene)
    chrom_b = setNames(b$chrom, b$interval_gene)
    expect_equal(chrom_a[a$interval_gene], chrom_b[a$interval_gene],
                 info = paste("chromosome disagreement in panel", nm))
    for (p in list(a, b)) {
      expect_false(any(is.na(p$start) | is.na(p$end)))
      expect_true(all(p$end >= p$start))
    }
  }
})

test_that("the sarcoma panel carries its fusion partners with hg38 coordinates", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  p = resolve_gene_panel("sarcoma", assets, "hg38")
  expect_true(p$has_coords)
  expect_equal(p$reference, "hg38")
  expect_equal(length(p$genes), 140L)
  # Canonical symbols, not the aliases the source list used.
  expect_true(all(c("EWSR1", "SS18", "MRTFB", "OGA", "KDR", "FLT4", "ERBB2",
                    "H3-3A", "POU2AF3", "DUX4L10") %in% p$genes))
  expect_false(any(c("MKL2", "MGEA5", "VEGFR2", "VEGFR3", "HER2", "SYT",
                     "H3F3A") %in% p$genes))
  # EWSR1's hg38 span, so a T2T file swapped in here would fail.
  expect_equal(p$chrom[p$interval_gene == "EWSR1"], "chr22")
  expect_equal(p$start[p$interval_gene == "EWSR1"], 29268009L)
})

test_that("load_all_gene_panels skips a panel that ships only for another reference", {
  repo_root = dirname(dirname(getwd()))
  # lymphoid ships for hg38 and t2t only, so an unrelated reference gets nothing:
  # offering it would mean matching another genome's coordinates.
  expect_message(panels <- load_all_gene_panels(file.path(repo_root, "assets"), "hs1"),
                 "ships only for reference")
  expect_false("lymphoid" %in% names(panels))
})

test_that("load_all_gene_panels panel names are lowercase filenames without extension", {
  repo_root = dirname(dirname(getwd()))
  panels = load_all_gene_panels(file.path(repo_root, "assets"), "hg38")
  expect_true(all(names(panels) == tolower(names(panels))))
})

test_that("resolve_gene_panel treats 'none' as no filtering", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  expect_null(resolve_gene_panel("none", assets))
  expect_null(resolve_gene_panel("NONE", assets))
  expect_null(resolve_gene_panel(" none ", assets))
  expect_null(resolve_gene_panel(NULL, assets))
  expect_null(resolve_gene_panel(NA_character_, assets))
})

test_that("resolve_gene_panel loads the builtin variant for the given reference", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  p = resolve_gene_panel("lymphoid", assets, "hg38")
  expect_true(p$has_coords)
  expect_equal(p$reference, "hg38")
  expect_true("MYC" %in% p$genes)
  # MYC's hg38 span, not the T2T one.
  expect_equal(p$start[p$interval_gene == "MYC"], 127735434L)
  expect_equal(resolve_gene_panel("lymphoid", assets, "t2t")$reference, "t2t")
})

test_that("resolve_gene_panel loads a symbol-only panel from a file path", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("gene", "BRCA1", "BRCA2"), tsv)
  on.exit(unlink(tsv))
  p = resolve_gene_panel(tsv, assets, "hg38")
  expect_false(p$has_coords)
  expect_equal(p$genes, c("BRCA1", "BRCA2"))
  expect_null(panel_intervals(p))
})

test_that("resolve_gene_panel errors on an unknown panel rather than silently unfiltering", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  expect_error(resolve_gene_panel("nonsense", assets, "hg38"), "Gene panel not found")
})

test_that("a headerless one-column panel keeps its first symbol", {
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("MYC", "BCL2"), tsv)
  on.exit(unlink(tsv))
  expect_setequal(load_gene_panel(tsv)$genes, c("MYC", "BCL2"))
})

test_that("a panel with a partial coordinate set is an error, not a silent downgrade", {
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("gene\tchrom\tstart", "MYC\tchr8\t127735434"), tsv)
  on.exit(unlink(tsv))
  expect_error(load_gene_panel(tsv), "but not all of")
})

test_that("a coordinate panel declaring the wrong reference is an error", {
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("# reference: hg38",
               "gene\tchrom\tstart\tend", "MYC\tchr8\t127735434\t127742951"), tsv)
  on.exit(unlink(tsv))
  expect_error(load_gene_panel(tsv, "t2t"), "declares reference")
  p = load_gene_panel(tsv, "hg38")
  expect_true(p$has_coords)
  expect_equal(p$reference, "hg38")
  # A reference column works the same way as the comment line.
  tsv2 = tempfile(fileext = ".tsv")
  writeLines(c("gene\tchrom\tstart\tend\treference",
               "MYC\tchr8\t127735434\t127742951\tgrch38"), tsv2)
  on.exit(unlink(tsv2), add = TRUE)
  expect_error(load_gene_panel(tsv2, "t2t"), "declares reference")
  expect_equal(load_gene_panel(tsv2, "hg38")$reference, "hg38")
})

test_that("a panel that declares no reference loads unverified rather than being guessed", {
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("gene\tchrom\tstart\tend", "MYC\tchr8\t127735434\t127742951"), tsv)
  on.exit(unlink(tsv))
  p = load_gene_panel(tsv, "t2t")
  expect_true(p$has_coords)
  expect_equal(p$reference, "")
})

test_that("a symbol-only panel is reference-agnostic", {
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("# reference: hg38", "gene", "MYC"), tsv)
  on.exit(unlink(tsv))
  expect_equal(load_gene_panel(tsv, "t2t")$genes, "MYC")
})

test_that("a coordinate panel with an unparseable coordinate is an error", {
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("gene\tchrom\tstart\tend", "MYC\tchr8\tnot_a_number\t127742951"), tsv)
  on.exit(unlink(tsv))
  expect_error(load_gene_panel(tsv), "missing or non-numeric")
})

test_that("js_col_index_map emits zero-based positions for the client-side filter", {
  expect_equal(js_col_index_map(c("id", "panel_hit")), '{"id":0,"panel_hit":1}')
  expect_equal(js_col_index_map(character(0)), "{}")
})
