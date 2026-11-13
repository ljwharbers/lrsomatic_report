# QC parsing and discovery. setup.R does not source these two files, so do it here.
source(file.path(repo_root, "R/parse_qc.R"))
source(file.path(repo_root, "R/locate_outputs.R"))

qc_fixture = function(name) file.path(repo_root, "tests", "fixtures", "qc", name)

# Build a sample tree from relative paths; each file gets `content` (default: one line)
make_tree = function(paths, content = "x") {
  root = tempfile("qc_tree_")
  for (p in paths) {
    f = file.path(root, p)
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeLines(content, f)
  }
  root
}
rel = function(path, root) if (is.null(path)) NULL else sub(paste0(root, "/"), "", path, fixed = TRUE)

# --- cramino -----------------------------------------------------------------

test_that("cramino 2.0: N50 is the full-read-length row, not 'N50 aligned'", {
  cr = parse_cramino(qc_fixture("cramino_v2.txt"))
  expect_equal(cr$n50, 25774)
  expect_equal(cr$n50_aligned, 25677)
  expect_equal(cr$n75, 21998)
  expect_equal(cr$n75_aligned, 21904)
  # "Yield [Gb]" must not be confused with "Yield [Gb] (>25kb)"
  expect_equal(cr$yield_gb, 76.27)
  expect_equal(cr$mapped_pct, 100)
  expect_identical(cr$n_reads, 2993079L)
  expect_equal(cr$median_identity, 99.67)
})

test_that("cramino 2.0: exact keys do not depend on row order", {
  # Put the aligned twins first: an unanchored grepl("N50") would now read 25677
  lines = readLines(qc_fixture("cramino_v2.txt"))
  aligned = grepl(" aligned\t", lines)
  f = tempfile(fileext = ".txt")
  writeLines(c(lines[aligned], lines[!aligned]), f)
  cr = parse_cramino(f)
  expect_equal(cr$n50, 25774)
  expect_equal(cr$n50_aligned, 25677)
  expect_equal(cr$yield_gb, 76.27)
})

test_that("cramino 1.x: no aligned rows, '% from total alignments' still read", {
  cr = parse_cramino(qc_fixture("cramino_v1.txt"))
  expect_equal(cr$n50, 25521)
  expect_true(is.na(cr$n50_aligned))
  expect_equal(cr$yield_gb, 76.48)
  expect_equal(cr$mapped_pct, 86.72)
  expect_identical(cr$n_reads, 2992478L)
})

test_that("cramino: missing or empty file gives every field as NA", {
  for (cr in list(parse_cramino(NULL), parse_cramino(tempfile()),
                  parse_cramino({f = tempfile(); file.create(f); f}))) {
    expect_setequal(names(cr), names(CRAMINO_KEYS))
    expect_true(all(vapply(cr, is.na, logical(1))))
  }
})

# --- flagstat ------------------------------------------------------------------

test_that("flagstat: 'mapped' and 'primary mapped' are not conflated", {
  fs = parse_flagstat(qc_fixture("flagstat_samtools113.txt"))
  expect_equal(fs$total, 3579257)
  expect_equal(fs$mapped, 3579197)
  expect_equal(fs$primary_mapped, 2993079)
  expect_equal(fs$primary, 2993139)
  expect_equal(fs$secondary, 552284)
  expect_equal(fs$supplementary, 33834)
  expect_equal(fs$duplicate, 0)
  expect_equal(fs$paired, 0)
})

test_that("flagstat: pre-1.13 output (no primary lines) still parses", {
  f = tempfile()
  writeLines(c("100 + 0 in total (QC-passed reads + QC-failed reads)",
               "5 + 0 secondary", "3 + 0 supplementary", "0 + 0 duplicates",
               "98 + 0 mapped (98.00% : N/A)", "0 + 0 paired in sequencing"), f)
  fs = parse_flagstat(f)
  expect_equal(fs$mapped, 98)
  expect_null(fs$primary_mapped)
  expect_equal(parse_flagstat(NULL), list())
})

# --- mosdepth / samtools stats robustness -------------------------------------------

test_that("0-byte stub QC files parse to empty results, not errors", {
  root = tempfile(); dir.create(root)
  f = file.path(root, "empty"); file.create(f)
  m = parse_mosdepth_summary(f)
  expect_true(is.na(m$mean_depth))
  expect_equal(nrow(m$table), 0)
  expect_null(parse_mosdepth_dist(f))
  expect_null(parse_samtools_stats(f))
})

# --- discovery -------------------------------------------------------------------

test_that("two same-named cramino files: cramino_aln wins over cramino_ubam", {
  root = make_tree(c("qc/tumor/cramino_aln/S_tumor_cramino.txt",
                     "qc/tumor/cramino_ubam_rep1/S_tumor_cramino.txt",
                     "qc/tumor/cramino_ubam_rep2/S_tumor_cramino.txt"))
  on.exit(unlink(root, recursive = TRUE))
  o = locate_outputs(root, "S")
  expect_equal(rel(o$cramino, root), "qc/tumor/cramino_aln/S_tumor_cramino.txt")
})

test_that("cramino: without a cramino_aln dir, a non-ubam hit beats the ubam one", {
  # "a_..." sorts before "cramino_ubam..." only by accident; the ranking must not depend on that
  root = make_tree(c("qc/tumor/cramino_ubam_rep1/S_tumor_cramino.txt",
                     "qc/tumor/z_other/S_tumor_cramino.txt"))
  on.exit(unlink(root, recursive = TRUE))
  expect_equal(rel(locate_outputs(root, "S")$cramino, root), "qc/tumor/z_other/S_tumor_cramino.txt")
})

test_that("published layout, current mosdepth naming, matched run", {
  # Mirrors a real stub run: qc/<type>/<tool>/, mosdepth prefix "<id>_<type>", no ".mosdepth" infix
  paths = c()
  for (t in c("tumor", "normal")) paths = c(paths,
    sprintf("qc/%s/cramino_aln/S_%s_cramino.txt", t, t),
    sprintf("qc/%s/cramino_ubam_rep1/S_%s_cramino.txt", t, t),
    sprintf("qc/%s/mosdepth/S_%s.summary.txt", t, t),
    sprintf("qc/%s/mosdepth/S_%s.global.dist.txt", t, t),
    sprintf("qc/%s/mosdepth/S_%s.region.dist.txt", t, t),
    sprintf("qc/%s/samtools/S_%s.stats", t, t),
    sprintf("qc/%s/samtools/S_%s.idxstats", t, t),
    sprintf("qc/%s/samtools/S_%s.flagstat", t, t))
  # A non-mosdepth ".summary.txt" outside qc/ must not be taken for the coverage summary
  paths = c("ecdna/coral/S.summary.txt", paths)
  root = make_tree(paths)
  on.exit(unlink(root, recursive = TRUE))
  o = locate_outputs(root, "S")
  expect_true(o$has_normal)
  expect_equal(o$mode, "matched")
  expect_equal(rel(o$mosdepth_summary, root), "qc/tumor/mosdepth/S_tumor.summary.txt")
  expect_equal(rel(o$mosdepth_dist, root), "qc/tumor/mosdepth/S_tumor.global.dist.txt")
  expect_equal(rel(o$cramino, root), "qc/tumor/cramino_aln/S_tumor_cramino.txt")
  expect_equal(rel(o$samtools_stats, root), "qc/tumor/samtools/S_tumor.stats")
  expect_equal(rel(o$flagstat, root), "qc/tumor/samtools/S_tumor.flagstat")
  expect_equal(rel(o$normal_mosdepth_summary, root), "qc/normal/mosdepth/S_normal.summary.txt")
  expect_equal(rel(o$normal_mosdepth_dist, root), "qc/normal/mosdepth/S_normal.global.dist.txt")
  expect_equal(rel(o$normal_cramino, root), "qc/normal/cramino_aln/S_normal_cramino.txt")
  expect_equal(rel(o$normal_samtools_stats, root), "qc/normal/samtools/S_normal.stats")
  expect_equal(rel(o$normal_flagstat, root), "qc/normal/samtools/S_normal.flagstat")
})

test_that("older mosdepth naming ('.mosdepth.' infix) is still found, tumour-only", {
  root = make_tree(c("qc/tumor/mosdepth/S_tumor.mosdepth.summary.txt",
                     "qc/tumor/mosdepth/S_tumor.mosdepth.global.dist.txt"))
  on.exit(unlink(root, recursive = TRUE))
  o = locate_outputs(root, "S")
  expect_equal(rel(o$mosdepth_summary, root), "qc/tumor/mosdepth/S_tumor.mosdepth.summary.txt")
  expect_equal(rel(o$mosdepth_dist, root), "qc/tumor/mosdepth/S_tumor.mosdepth.global.dist.txt")
  expect_false(o$has_normal)
  expect_equal(o$mode, "tumour-only")
})

test_that("pipeline-staged flat layout (qc/<type>/ symlinks, no tool subdirs)", {
  src = make_tree(c("S_tumor.summary.txt", "S_tumor.global.dist.txt", "S_tumor_cramino.txt",
                    "S_tumor.stats", "S_tumor.flagstat", "S_normal.summary.txt",
                    "S_normal_cramino.txt"))
  root = tempfile("staged_")
  dir.create(file.path(root, "qc", "tumor"), recursive = TRUE)
  dir.create(file.path(root, "qc", "normal"), recursive = TRUE)
  for (f in list.files(src)) {
    side = if (grepl("_normal", f)) "normal" else "tumor"
    file.symlink(file.path(src, f), file.path(root, "qc", side, f))
  }
  on.exit(unlink(c(root, src), recursive = TRUE))
  o = locate_outputs(root, "S")
  expect_equal(rel(o$mosdepth_summary, root), "qc/tumor/S_tumor.summary.txt")
  expect_equal(rel(o$mosdepth_dist, root), "qc/tumor/S_tumor.global.dist.txt")
  expect_equal(rel(o$cramino, root), "qc/tumor/S_tumor_cramino.txt")
  expect_equal(rel(o$samtools_stats, root), "qc/tumor/S_tumor.stats")
  expect_equal(rel(o$flagstat, root), "qc/tumor/S_tumor.flagstat")
  expect_equal(rel(o$normal_cramino, root), "qc/normal/S_normal_cramino.txt")
  expect_equal(rel(o$normal_mosdepth_summary, root), "qc/normal/S_normal.summary.txt")
  expect_null(o$normal_flagstat)
  expect_true(o$has_normal)
})

test_that("a sample_dir whose own path contains /normal/ is not read as normal-side", {
  base = tempfile("qc_")
  root = file.path(base, "normal", "S")
  dir.create(file.path(root, "qc", "tumor"), recursive = TRUE)
  writeLines("x", file.path(root, "qc", "tumor", "S_tumor_cramino.txt"))
  on.exit(unlink(base, recursive = TRUE))
  o = locate_outputs(root, "S")
  expect_false(is.null(o$cramino))
  expect_false(o$has_normal)
})
