test_that("whatshap_totals reads the ALL row when whatshap wrote one", {
  w = list(
    per_chrom = data.table(chromosome = c("chr1", "chr2"),
                           variants = c(10, 20), phased = c(5, 10)),
    all = list(chromosome = "ALL", variants = 30, phased = 15, phased_fraction = 0.5)
  )
  got = whatshap_totals(w)
  expect_equal(got$phased, 15)
  expect_equal(got$variants, 30)
  expect_equal(got$fraction, 0.5)
})

test_that("whatshap_totals sums the per-chromosome rows when there is no ALL row", {
  # Silent N/A before: the Phasing table renders fine while the header card has nothing
  w = list(
    per_chrom = data.table(chromosome = c("chr1", "chr2"),
                           variants = c(10, 30), phased = c(5, 10),
                           phased_fraction = c(0.5, 1 / 3)),
    all = NULL
  )
  got = whatshap_totals(w)
  expect_equal(got$phased, 15)
  expect_equal(got$variants, 40)
  # Recomputed from the totals, not averaged over the per-chromosome fractions
  expect_equal(got$fraction, 15 / 40)
})

test_that("whatshap_totals recomputes the fraction when the column is absent", {
  # This shape used to error in the card: is.na(NULL) is logical(0), not FALSE
  w = list(per_chrom = data.table(chromosome = "chr1"),
           all = list(chromosome = "ALL", variants = 200, phased = 150))
  got = whatshap_totals(w)
  expect_equal(got$fraction, 0.75)
})

test_that("whatshap_totals returns NULL when the summed columns are entirely missing", {
  # whatshap leaves the field empty for a chromosome with no het variants. na.rm = TRUE made
  # the sum 0, which passed the NULL guard and put a confident "Phased variants 0" in the
  # header over a file that reported nothing at all.
  w = list(
    per_chrom = data.table(chromosome = c("chr1", "chr2"),
                           variants = c("", ""), phased = c("", "")),
    all = NULL
  )
  expect_null(whatshap_totals(w))
})

test_that("whatshap_totals returns NULL when there is nothing to total", {
  expect_null(whatshap_totals(NULL))
  expect_null(whatshap_totals(list(per_chrom = data.table(), all = NULL)))
  # An ALL row carrying neither a count nor a fraction is not a number to show
  expect_null(whatshap_totals(list(per_chrom = data.table(chromosome = "chr1"),
                                   all = list(chromosome = "ALL", blocks = 4))))
})

test_that("the whatshap section parses a stats TSV with a #-prefixed header", {
  f = tempfile(fileext = ".tsv")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(
    "#sample\tchromosome\tfile_name\tvariants\tphased\tphased_fraction",
    "S1\tchr1\tgermline_smallvariants.vcf.gz\t100\t80\t0.8",
    "S1\tALL\tgermline_smallvariants.vcf.gz\t100\t80\t0.8"
  ), f)

  parsed = SECTIONS$whatshap$parse(list(stats_tsv = f), list())
  expect_equal(names(parsed$per_chrom)[1:2], c("sample", "chromosome"))
  expect_equal(nrow(parsed$per_chrom), 1L)
  expect_equal(parsed$vcf, "germline_smallvariants.vcf.gz")
  expect_equal(whatshap_totals(parsed)$phased, 80)
})

test_that("the whatshap section returns NULL when the stats file is absent", {
  expect_null(SECTIONS$whatshap$parse(list(stats_tsv = NULL), list()))
  expect_null(SECTIONS$whatshap$parse(list(stats_tsv = "/no/such/file.tsv"), list()))
})

test_that("the whatshap section returns NULL on a 0-byte stats file (stub runs)", {
  f = tempfile(fileext = ".tsv"); file.create(f)
  expect_null(suppressMessages(SECTIONS$whatshap$parse(list(stats_tsv = f), list())))
})

test_that("whatshap stats are found in the published qc/whatshap_stats/ layout", {
  root = tempfile("ws_")
  dir.create(file.path(root, "qc", "whatshap_stats"), recursive = TRUE)
  for (ext in c("tsv", "gtf", "log")) writeLines("x", file.path(root, "qc", "whatshap_stats", paste0("S_whatshap_stats.", ext)))
  on.exit(unlink(root, recursive = TRUE))
  got = SECTIONS$whatshap$locate(root, "S")$stats_tsv
  expect_equal(basename(got), "S_whatshap_stats.tsv")
})

test_that("the phased fraction is of heterozygous variants, as whatshap reports it (REP-4)", {
  # A2 CCS15-ONT's ALL row: 92.6 % = phased / heterozygous, not phased / variants (56.6 %)
  w = list(per_chrom = data.table(chromosome = "chr1"),
           all = list(chromosome = "ALL", variants = 5063778, phased = 2864966,
                      heterozygous_variants = 3092264, phased_fraction = 0.926494633058497))
  got = whatshap_totals(w)
  expect_equal(got$heterozygous, 3092264)
  expect_equal(got$denominator, "heterozygous")
  expect_equal(got$fraction, 0.926494633058497)
  # Summed rows recompute against the heterozygous count too
  s = whatshap_totals(list(all = NULL, per_chrom = data.table(
    chromosome = c("chr1", "chr2"), variants = c(100, 100), phased = c(40, 50),
    heterozygous_variants = c(50, 50))))
  expect_equal(s$fraction, 0.9)
  # Without a heterozygous count the fraction falls back to variants, and says so
  v = whatshap_totals(list(all = list(chromosome = "ALL", variants = 200, phased = 150), per_chrom = NULL))
  expect_equal(v$denominator, "variants")
  expect_equal(v$fraction, 0.75)
})
