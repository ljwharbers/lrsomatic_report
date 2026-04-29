# Minimal synthetic VEP data (mimics parse_vep_text output)
make_vep = function() {
  data.table(
    chrom        = c("chr1", "chr1", "chr2"),
    pos          = c(100L, 200L, 300L),
    ref          = c("A", "C", "G"),
    alt          = c("T", "G", "A"),
    symbol       = c("MYC", "TP53", "KRAS"),
    gene_id      = c("ENSG001", "ENSG002", "ENSG003"),
    consequence  = c("missense_variant", "stop_gained", "synonymous_variant"),
    impact       = c("MODERATE", "HIGH", "LOW"),
    hgvsp        = c("p.A1T", "p.Q2*", NA_character_),
    existing     = c(NA_character_, NA_character_, NA_character_),
    dbsnp        = c(NA_character_, NA_character_, NA_character_),
    cosmic       = c(NA_character_, NA_character_, NA_character_),
    sift         = c(NA_character_, NA_character_, NA_character_),
    polyphen     = c(NA_character_, NA_character_, NA_character_)
  )
}

test_that("build_variant_table with NULL panel returns all variants", {
  vep = make_vep()
  result = build_variant_table(vep, list(), gene_panel = NULL)
  expect_equal(nrow(result), 3L)
  expect_true(all(c("MYC", "TP53", "KRAS") %in% result$symbol))
})

test_that("build_variant_table with panel filters to panel genes", {
  vep = make_vep()
  result = build_variant_table(vep, list(), gene_panel = c("MYC", "TP53"))
  expect_equal(nrow(result), 2L)
  expect_false("KRAS" %in% result$symbol)
})

test_that("build_variant_table with empty panel returns empty table", {
  vep = make_vep()
  result = build_variant_table(vep, list(), gene_panel = character(0))
  expect_equal(nrow(result), 0L)
})
