make_sv_tsv = function() {
  data.table(
    ID           = c("SV1", "SV2", "SV3"),
    SVTYPE       = c("DEL", "DUP", "INV"),
    DETAILED_TYPE = c("DEL", "DUP", "INV"),
    START_CHROM  = c("chr1", "chr2", "chr3"),
    START_POS    = c(1000L, 2000L, 3000L),
    END_CHROM    = c("chr1", "chr2", "chr3"),
    END_POS      = c(5000L, 6000L, 7000L),
    SV_LEN       = c(4000L, 4000L, 4000L),
    VAF          = c(0.4, 0.5, 0.3),
    NHL_GENE_HITS = c("MYC;BCL2", "TP53", NA_character_),
    gene_hits    = c("MYC;BCL2", "TP53", NA_character_)
  )
}

test_that("build_sv_table with NULL panel returns one row per SV", {
  sv = make_sv_tsv()
  result = build_sv_table(sv, gene_panel = NULL)
  expect_equal(nrow(result), 3L)
})

test_that("build_sv_table with panel filters and explodes rows", {
  sv = make_sv_tsv()
  result = build_sv_table(sv, gene_panel = c("MYC", "TP53"))
  # SV1 contributes 1 row (MYC matches; BCL2 does not), SV2 contributes 1 row (TP53)
  expect_equal(nrow(result), 2L)
  expect_true(all(result$gene %in% c("MYC", "TP53")))
})

test_that("build_sv_table with empty panel returns empty table", {
  sv = make_sv_tsv()
  result = build_sv_table(sv, gene_panel = character(0))
  expect_equal(nrow(result), 0L)
})
