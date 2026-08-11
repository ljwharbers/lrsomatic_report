test_that("setup loads without error", {
  expect_true(is.function(load_gene_panel))
})

test_that("load_all_gene_panels returns named list of character vectors", {
  # Use the real assets dir
  # When testthat runs, getwd() is in tests/testthat, so we need to go up two levels to repo root
  repo_root = dirname(dirname(getwd()))
  panels = load_all_gene_panels(file.path(repo_root, "assets"))
  expect_type(panels, "list")
  expect_true(length(panels) >= 1)
  expect_true("lymphoid" %in% names(panels))
  expect_type(panels[["lymphoid"]], "character")
  expect_true(length(panels[["lymphoid"]]) > 0)
})

test_that("load_all_gene_panels panel names are lowercase filenames without extension", {
  repo_root = dirname(dirname(getwd()))
  panels = load_all_gene_panels(file.path(repo_root, "assets"))
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

test_that("resolve_gene_panel loads a builtin panel by name", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  genes = resolve_gene_panel("lymphoid", assets)
  expect_type(genes, "character")
  expect_true(length(genes) > 0)
  expect_true("MYC" %in% genes)
})

test_that("resolve_gene_panel loads a panel from a file path", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  tsv = tempfile(fileext = ".tsv")
  writeLines(c("gene", "BRCA1", "BRCA2"), tsv)
  on.exit(unlink(tsv))
  expect_equal(resolve_gene_panel(tsv, assets), c("BRCA1", "BRCA2"))
})

test_that("resolve_gene_panel errors on an unknown panel rather than silently unfiltering", {
  assets = file.path(dirname(dirname(getwd())), "assets")
  expect_error(resolve_gene_panel("nonsense", assets), "Gene panel not found")
})

test_that("filter_by_gene_panel returns dt unchanged when panel_genes is NULL", {
  dt = data.table(gene = c("MYC", "FOO"), n = 1:2)
  expect_identical(filter_by_gene_panel(dt, NULL), dt)
})

test_that("filter_by_gene_panel filters to the panel, and an empty panel keeps nothing", {
  dt = data.table(gene = c("MYC", "FOO"), n = 1:2)
  expect_equal(filter_by_gene_panel(dt, "MYC")$gene, "MYC")
  expect_equal(nrow(filter_by_gene_panel(dt, character(0))), 0L)
})
