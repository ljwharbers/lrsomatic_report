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
  result = build_variant_table(vep, NULL, gene_panel = NULL)
  expect_equal(nrow(result), 3L)
  expect_true(all(c("MYC", "TP53", "KRAS") %in% result$symbol))
})

test_that("build_variant_table with panel filters to panel genes", {
  vep = make_vep()
  result = build_variant_table(vep, NULL, gene_panel = c("MYC", "TP53"))
  expect_equal(nrow(result), 2L)
  expect_false("KRAS" %in% result$symbol)
})

test_that("build_variant_table with empty panel returns empty table", {
  vep = make_vep()
  result = build_variant_table(vep, NULL, gene_panel = character(0))
  expect_equal(nrow(result), 0L)
})

# --- vaf_provenance() -------------------------------------------------------

# Write a minimal gzipped VCF header (plus one record, so the scan has to stop on its own)
make_vcf_gz = function(path, source_line = NULL) {
  lines = c("##fileformat=VCFv4.2")
  if (!is.null(source_line)) lines = c(lines, paste0("##source=", source_line))
  lines = c(lines,
            "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
            "chr1\t100\t.\tA\tT\t30\tPASS\t.\tGT:AF\t0/1:0.4")
  con = gzfile(path, "wb")
  writeLines(lines, con)
  close(con)
  path
}

test_that("vaf_provenance returns NULL when there is nothing to describe", {
  expect_null(vaf_provenance(NULL))
  expect_null(vaf_provenance(character(0)))
  expect_null(vaf_provenance(NA_character_))
  expect_null(vaf_provenance(""))
})

test_that("vaf_provenance reports the path and the declared ##source", {
  d = withr::local_tempdir()
  f = make_vcf_gz(file.path(d, "somatic.vcf.gz"), source_line = "Clair3")

  prov = vaf_provenance(f, sample_dir = d)
  expect_equal(prov$paths, "somatic.vcf.gz")   # relative to sample_dir
  expect_equal(prov$sources, "Clair3")
})

test_that("vaf_provenance keeps the absolute path when it is outside sample_dir", {
  d = withr::local_tempdir()
  f = make_vcf_gz(file.path(d, "somatic.vcf.gz"), source_line = "Clair3")

  prov = vaf_provenance(f, sample_dir = file.path(d, "elsewhere"))
  expect_equal(prov$paths, f)
})

test_that("vaf_provenance handles several files and a missing ##source", {
  d = withr::local_tempdir()
  snvs  = make_vcf_gz(file.path(d, "snvs.vcf.gz"),  source_line = "ClairS")
  indel = make_vcf_gz(file.path(d, "indel.vcf.gz"), source_line = NULL)

  prov = vaf_provenance(c(snvs, indel), sample_dir = d)
  expect_equal(prov$paths, c("snvs.vcf.gz", "indel.vcf.gz"))
  expect_equal(prov$sources, "ClairS")   # deduped; the second file declares none
})

test_that("vaf_provenance does not fail on a path that does not exist", {
  prov = vaf_provenance("/nonexistent/somatic.vcf.gz")
  expect_equal(prov$paths, "/nonexistent/somatic.vcf.gz")
  expect_length(prov$sources, 0L)
})
