# Padfoot annotation (R/parse_padfoot.R). No real Padfoot output existed when this was written; the
# fixture follows the 26-column header of the pipeline's Padfoot module and its nf-test, keyed the
# way Padfoot keys SAVANA (the shared ID stem, ID_1 for ID_1_1/ID_1_2).
pf_fixture_dir = file.path(dirname(dirname(getwd())), "tests", "fixtures", "sv")
pf_savana_tsv  = file.path(pf_fixture_dir, "padfoot_savana_annotated_svs.tsv")
pf_savana_vcf  = file.path(pf_fixture_dir, "sample.classified.somatic.vcf")

write_padfoot = function(rows) {
  lines = readLines(pf_savana_tsv)
  path = tempfile(fileext = ".tsv")
  writeLines(c(lines[1], rows), path)
  path
}

test_that("parse_padfoot keeps pos_1 and Pos_1 apart and merges both sides", {
  pf = parse_padfoot(pf_savana_tsv)
  expect_equal(names(pf), c("sv_id", "padfoot_impact", "padfoot_score", "padfoot_role",
                            "padfoot_genes", "padfoot_type", "chrom_1", "pos_1",
                            "chrom_2", "pos_2"))
  r = pf[sv_id == "ID_25860"]
  # pos_1 is the coordinate, Pos_1 the position in the gene
  expect_equal(r$pos_1, 635796L)
  expect_equal(r$padfoot_genes, "GENEA (intron3), GENEA (exon10)")
  expect_equal(r$padfoot_role, "TSG, fusion")
  expect_equal(r$padfoot_impact, "HIGH_ONCO")
  expect_equal(r$padfoot_score, 8.5)
  expect_equal(r$padfoot_type, "between_exons")
  # A side with no gene contributes nothing; "NA" cells are missing values
  expect_equal(pf[sv_id == "ID_27202", padfoot_genes], "GENEB (upstream)")
  expect_equal(pf[sv_id == "ID_26555", padfoot_role], NA_character_)
  expect_true(is.na(pf[sv_id == "ID_99999", padfoot_impact]))
})

test_that("a missing Padfoot file is NULL and an empty one has no rows", {
  expect_null(parse_padfoot(NULL))
  expect_null(parse_padfoot(tempfile()))
  e = tempfile(); file.create(e); on.exit(unlink(e))
  expect_equal(nrow(parse_padfoot(e)), 0L)
})

test_that("SAVANA rows join on the stripped id stem, chosen once for the table", {
  t = build_sv_table_from_vep(parse_sv_vcf(pf_savana_vcf, "savana"), NULL)
  j = join_padfoot(t, parse_padfoot(pf_savana_tsv))
  expect_equal(j$key, "id stem")
  expect_equal(j$n_joined, 3L)
  expect_equal(j$n_rows, 12L)
  expect_equal(j$table[id == "ID_25860_1", padfoot_impact], "HIGH_ONCO")
  expect_equal(j$table[id == "ID_26555_1", padfoot_genes], "KAZN (intron3)")
  expect_true(is.na(j$table[id == "ID_29721_1", padfoot_impact]))
  expect_true(all(PADFOOT_COLS %in% names(j$table)))
})

severus_pair_table = function() data.table(
  id = c("severus_BND1_1", "severus_DEL1"), id_b = c("severus_BND1_2", NA),
  svclass = c("translocation", "DEL"), svtype = c("BND", "DEL"),
  chrom_a = c("chr1", "chr2"), pos_a = c(1000L, 5000L),
  chrom_b = c("chr9", "chr2"), pos_b = c(7000L, 9000L))

pf_row = function(id, c1, p1, c2, p2, impact)
  paste(c(id, c1, p1, c2, p2, "5", "0.3", "+-", "1", impact, "NA", "MYC", "exon2", "+", "NA",
          "NA", "NA", "NA", "NA", "NA", "exonic", "0", "False", "False", "NA", "NA"),
        collapse = "\t")

test_that("Severus rows join on their own id (Padfoot reports a pair under its _1 record)", {
  f = write_padfoot(c(pf_row("severus_BND1_1", "chr1", 1000, "chr9", 7000, "HIGH_TSG"),
                      pf_row("severus_DEL1", "chr2", 5000, "chr2", 9000, "LOW")))
  on.exit(unlink(f))
  j = join_padfoot(severus_pair_table(), parse_padfoot(f))
  expect_equal(j$key, "id")
  expect_equal(j$n_joined, 2L)
})

test_that("the mate id is tried when Padfoot kept the other record", {
  f = write_padfoot(pf_row("severus_BND1_2", "chr9", 7000, "chr1", 1000, "HIGH_ONCO"))
  on.exit(unlink(f))
  j = join_padfoot(severus_pair_table(), parse_padfoot(f))
  expect_equal(j$key, "mate id")
  expect_equal(j$table[id == "severus_BND1_1", padfoot_impact], "HIGH_ONCO")
})

test_that("no id agreeing falls back to the loci, in either orientation", {
  f = write_padfoot(c(pf_row("other_1", "9", 7000, "1", 1000, "HIGH_ONCO"),   # reversed, bare contig
                      pf_row("other_2", "chr2", 5000, "chr2", 9000, "LOW")))
  on.exit(unlink(f))
  j = join_padfoot(severus_pair_table(), parse_padfoot(f))
  expect_equal(j$key, "locus")
  expect_equal(j$n_joined, 2L)
  expect_equal(j$table[id == "severus_BND1_1", padfoot_impact], "HIGH_ONCO")
})

test_that("nothing joining leaves the columns empty and says so", {
  f = write_padfoot(pf_row("x", "chr5", 1, "chr5", 2, "LOW"))
  on.exit(unlink(f))
  j = join_padfoot(severus_pair_table(), parse_padfoot(f))
  expect_true(is.na(j$key))
  expect_equal(j$n_joined, 0L)
  expect_true(all(is.na(j$table$padfoot_impact)))
})

test_that("the SV section adds Padfoot columns only when a Padfoot file exists", {
  d = tempfile("sample_")
  dir.create(file.path(d, "variants/savana"), recursive = TRUE)
  dir.create(file.path(d, "padfoot/savana"), recursive = TRUE)
  file.copy(pf_savana_vcf, file.path(d, "variants/savana/sample.classified.somatic.vcf"))
  on.exit(unlink(d, recursive = TRUE))

  res = SECTIONS$sv$parse(SECTIONS$sv$locate(d, "sample"), list())
  expect_false(any(PADFOOT_COLS %in% names(res$table)))

  file.copy(pf_savana_tsv, file.path(d, "padfoot/savana/annotated_svs.tsv"))
  inp = SECTIONS$sv$locate(d, "sample")
  expect_equal(inp$callers$savana$padfoot$combo, "savana")
  expect_null(inp$callers$severus$padfoot)
  res = SECTIONS$sv$parse(inp, list())
  expect_true(res$has_padfoot)
  expect_true(all(PADFOOT_COLS %in% names(res$table)))
  expect_equal(res$padfoot$savana$key, "id stem")
  expect_equal(res$padfoot$savana$n_joined, 3L)
  expect_equal(res$padfoot$savana$n_padfoot, 4L)
})
