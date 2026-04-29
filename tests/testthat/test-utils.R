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
