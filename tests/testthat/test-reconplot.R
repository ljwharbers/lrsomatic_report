# ReConPlot discovery (R/sections/reconplot.R). PNG content is irrelevant here: empty files stand
# in, and the pipeline's staging (individual symlinked files in a mirrored tree) is exercised.

make_reconplot_tree = function(sample_id = "S1", symlink = FALSE) {
  d   = tempfile("rp_")
  src = tempfile("rp_src_")
  dir.create(src)
  add = function(rel) {
    dest = file.path(d, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (symlink) {
      s = file.path(src, gsub("/", "__", rel))
      file.create(s)
      file.symlink(s, dest)
    } else file.create(dest)
  }
  for (cb in c("savana", "severus_ascat")) {
    add(sprintf("reconplot/%s/genome_wide/%s_genome_wide.png", cb, sample_id))
    add(sprintf("reconplot/%s/genome_wide/%s_genome_wide.pdf", cb, sample_id))
    for (ch in c("chr10", "chr2", "chrX", "chr1"))
      add(sprintf("reconplot/%s/per_chromosome/%s_%s.png", cb, sample_id, ch))
    add(sprintf("reconplot/%s/focus/%s_MYC.png", cb, sample_id))
  }
  list(dir = d, src = src)
}

rp_ctx = list(chromosomes = paste0("chr", c(1:22, "X", "Y")), embed = "slim")

test_that("ReConPlot figures are found per combination, PDFs and focus plots ignored", {
  t = make_reconplot_tree()
  on.exit(unlink(c(t$dir, t$src), recursive = TRUE))
  s = SECTIONS$reconplot
  res = s$parse(s$locate(t$dir, "S1"), list(), rp_ctx)
  # Known combinations in their fixed order, not alphabetical
  expect_equal(names(res$combos), c("severus_ascat", "savana"))
  x = res$combos$severus_ascat
  expect_match(x$genome, "severus_ascat/genome_wide/S1_genome_wide\\.png$")
  # Display order from ctx$chromosomes
  expect_equal(names(x$chroms), c("chr1", "chr2", "chr10", "chrX"))
  expect_false(any(grepl("focus|\\.pdf$", c(x$genome, x$chroms))))
})

test_that("staged symlinked files are discovered", {
  skip_on_os("windows")
  t = make_reconplot_tree(symlink = TRUE)
  on.exit(unlink(c(t$dir, t$src), recursive = TRUE))
  s = SECTIONS$reconplot
  res = s$parse(s$locate(t$dir, "S1"), list(), rp_ctx)
  expect_equal(length(res$combos), 2L)
  expect_equal(length(res$combos$savana$chroms), 4L)
})

test_that("no ReConPlot output parses to NULL", {
  d = tempfile("rp_"); dir.create(file.path(d, "variants"), recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))
  s = SECTIONS$reconplot
  expect_null(s$parse(s$locate(d, "S1"), list(), rp_ctx))
})
