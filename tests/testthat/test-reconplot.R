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

# ---- which per-chromosome figures slim mode embeds -----------------------------------------

test_that("ReConPlot chromosomes: ShatterSeek and Padfoot HIGH always, findings per list", {
  ss = data.table(chrom = c("7", "chr3"), confidence = c("high", "none"))
  # Plain HIGH (chr20) is too common to select by; HIGH_ONCO is the cancer-gene rating
  # Both mates of the translocation, in Padfoot's "GENE (where)" notation
  sv = data.table(chrom_a = c("chr22", "chr12", "chr1", "chr20"),
                  chrom_b = c("chr12", "chr22", "chr1", "chr20"),
                  padfoot_impact = c("HIGH_ONCO", "HIGH_ONCO", "LOW", "HIGH"),
                  padfoot_genes = c("EWSR1 (intron7), ATF1 (intron3)",
                                    "ATF1 (intron3), EWSR1 (intron7)", NA, "MACROD2"))
  f = data.table(gene = c("EWSR1–ATF1", "CDKN2A", "CDKN2A", "KMT2D", "MYC"),
                 genes = c("EWSR1,ATF1", "CDKN2A", "CDKN2A", "KMT2D", "MYC"),
                 type = c("SV", "SV", "CNA", "small variant", "SV"),
                 tier = c(1L, 2L, 1L, 1L, 3L),
                 panels = c("sarcoma", "sarcoma,lymphoid", "lymphoid", "lymphoid", "lymphoid"),
                 alteration = c("EWSR1–ATF1 translocation chr12–chr22 (SAVANA, Severus)",
                                "CDKN2A deletion chr9:21.97–24.84 Mb (SAVANA)",
                                "homozygous deletion (CN 0, partial gene)", "p.R1*", "near"),
                 locus = c("chr12–chr22", "chr9:21.90–22.01 Mb", "chr9:21,967,752–21,995,301",
                           "chr12:49,020,000", "chr8:127.7 Mb"))
  r = reconplot_chromosomes(f, sv, ss)
  setkey(r, chrom)
  # chr3's call has confidence none, chr1's SV is not HIGH; small variants and tier 3 never select
  expect_setequal(r$chrom, c("chr7", "chr22", "chr12", "chr9"))
  expect_true(r["chr7"]$always)
  expect_true(r["chr22"]$always)
  # One Padfoot reason per rating, symbols only, mates collapsed; alterations cut at the locus
  expect_equal(r["chr22"]$reasons,
               "Padfoot HIGH_ONCO: EWSR1, ATF1 · EWSR1–ATF1 translocation (tier 1)")
  expect_equal(r["chr9"]$reasons,
               "CDKN2A deletion (tier 2) · CDKN2A homozygous deletion (tier 1)")
  # chr9 is list-dependent; its lists are the union of both findings
  expect_false(r["chr9"]$always)
  expect_setequal(strsplit(r["chr9"]$panels, ",")[[1]], c("sarcoma", "lymphoid"))
  expect_equal(r["chr9"]$genes, "CDKN2A")
  # An always-shown chromosome carries no list (nothing to hide it by)
  expect_equal(r["chr22"]$panels, "")
  expect_equal(nrow(reconplot_chromosomes()), 0)
})

test_that("ReConPlot reasons are capped", {
  f = data.table(gene = paste0("G", 1:5), genes = paste0("G", 1:5), type = "CNA", tier = 1L,
                 panels = "a", alteration = "amplification", locus = "chr2:1–2")
  r = reconplot_chromosomes(f)
  expect_match(r$reasons, "2 more$")
})
