# Copy number (R/sections/cna.R + parsers in R/parse_ascat.R). Fixtures in tests/fixtures/cna/
# are hand-made in the formats Wakhan 0.5, ASCAT and SAVANA write (checked against pilot15).
cna_fx = file.path(repo_root, "tests", "fixtures", "cna")
cna_fx_path = function(...) file.path(cna_fx, ...)

# A panel as load_gene_panel() builds it (only the fields cna.R reads)
fake_panel = function(genes, chrom = character(0), start = integer(0), end = integer(0),
                      coord_genes = character(0), role_gene = character(0), role = character(0)) {
  list(name = "p", genes = genes, has_coords = length(chrom) > 0,
       all_chrom = chrom, all_start = as.integer(start), all_end = as.integer(end),
       all_gene = coord_genes, role_gene = role_gene, role = role)
}

seg_dt = function(chrom, start, end, a, b) {
  dt = data.table(chrom = chrom, start = as.numeric(start), end = as.numeric(end),
                  cn_a = as.numeric(a), cn_b = as.numeric(b))
  dt[, `:=`(total = cn_a + cn_b, minor = pmin(cn_a, cn_b), conf = NA_real_,
            sv_ids = NA_character_, conf_a = NA_real_, conf_b = NA_real_)]
  dt
}

# A fresh directory removed when the calling test ends: test-foundation.R lists tempdir()
# recursively and expects to find no CN output there (withr ships with testthat)
cna_tempdir = function(name, envir = parent.frame(2)) {
  root = tempfile("cna")
  d = file.path(root, name)
  dir.create(d, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE), envir = envir)
  d
}

# Lay the fixtures out as a sample directory: `mode` "copy" (real output dir, solution_rank_1
# a symlinked DIRECTORY onto the real solution dir, as Wakhan writes it) or "staged" (the
# pipeline's report task: every file an individual symlink, no directories linked)
make_sample_dir = function(mode = c("copy", "staged"), wakhan = TRUE, ascat = TRUE, savana = TRUE) {
  mode = match.arg(mode)
  d = cna_tempdir("S1")
  put = function(src, dst, link = FALSE) {
    dir.create(dirname(file.path(d, dst)), recursive = TRUE, showWarnings = FALSE)
    if (link) file.symlink(src, file.path(d, dst)) else file.copy(src, file.path(d, dst))
  }
  if (wakhan) {
    if (mode == "copy") {
      real = "wakhan/solution_2.02_0.57_0.92"
      put(cna_fx_path("wakhan/solution_rank_1/integer_profile.bed"), file.path(real, "integer_profile.bed"))
      put(cna_fx_path("wakhan/solution_rank_1/genes/genes_copynumber_states.bed"),
          file.path(real, "genes/genes_copynumber_states.bed"))
      file.symlink("solution_2.02_0.57_0.92", file.path(d, "wakhan", "solution_rank_1"))
      put(cna_fx_path("wakhan/solutions_ranks.tsv"), "wakhan/solutions_ranks.tsv")
    } else {
      for (f in c("solutions_ranks.tsv", "solution_rank_1/integer_profile.bed",
                  "solution_rank_1/genes/genes_copynumber_states.bed"))
        put(cna_fx_path("wakhan", f), file.path("wakhan", f), link = TRUE)
    }
  }
  if (ascat) for (f in c("S1.segments_raw.txt", "S1.purityploidy.txt", "S1.metrics.txt"))
    put(cna_fx_path("ascat", f), if (mode == "copy") file.path("ascat", f) else f, link = mode == "staged")
  if (savana) put(cna_fx_path("S1_fitted_purity_ploidy.tsv"),
                  "variants/savana/S1_fitted_purity_ploidy.tsv", link = mode == "staged")
  d
}

cna_ctx = function(sex = "female", panels = list()) {
  list(sample_id = "S1", reference = "t2t", sex = sex, all_panels = panels, embed = "slim",
       chromosomes = chromosomes_for_sex(sex),
       chrom_lens = setNames(rep(5e6, 24), c(paste0("chr", 1:22), "chrX", "chrY")))
}

test_that("cna_call applies the COSMIC thresholds on both sides of ploidy 2.7", {
  # ploidy 2: amp >= 5, homdel == 0, else vs expected 2
  low = cna_call(total = c(0, 1, 2, 3, 4, 5, 12), minor = 1, ploidy = 2, expected = 2)
  expect_equal(low$call, c("homozygous deletion", "loss", "neutral", "gain", "gain",
                           "amplification", "amplification"))
  # exactly at the split still counts as the low side
  expect_equal(cna_call(5, 1, 2.7, 3)$call, "amplification")
  # ploidy 4: amp >= 9, deletion < 4 - 2.7 = 1.3: CN 0 is homozygous, CN 1 deep (REP-7)
  high = cna_call(total = c(0, 1, 2, 3, 4, 5, 8, 9), minor = 1, ploidy = 4, expected = 4)
  expect_equal(high$call, c("homozygous deletion", "deep deletion", "loss", "loss",
                            "neutral", "gain", "gain", "amplification"))
  expect_true(is.na(cna_call(NA, 1, 2, 2)$call))
  expect_true(is.na(cna_call(2, 1, NA, 2)$call))
})

test_that("LOH needs a missing minor allele on a non-zero segment where two are expected", {
  r = cna_call(total = c(2, 0, 1, 1), minor = c(0, 0, 0, 0), ploidy = 2, expected = c(2, 2, 2, 1))
  expect_equal(r$loh, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(r$call, c("neutral", "homozygous deletion", "loss", "neutral"))
})

test_that("expected CN is halved on chrX/chrY in a male only", {
  expect_equal(cna_expected_cn(c("chr1", "chrX", "chrY"), 2.02, "male"), c(2, 1, 1))
  expect_equal(cna_expected_cn(c("chr1", "chrX"), 2.02, "female"), c(2, 2))
  expect_equal(cna_expected_cn(c("chr1", "chrX"), 3.94, "male"), c(4, 2))
  # One X copy in a male is neutral, in a female a loss
  segs = seg_dt("chrX", 1, 1e6, 1, 0)
  p = list(x = fake_panel("GX", "chrX", 100, 200, "GX"))
  m = cna_gene_calls(segs, p, 2, "male", "ascat")$genes
  f = cna_gene_calls(segs, p, 2, "female", "ascat")$genes
  expect_equal(m$call, "neutral"); expect_false(m$loh); expect_equal(m$expected, 1)
  expect_equal(f$call, "loss");    expect_true(f$loh)
})

test_that("a gene spanning several segments is judged on its extremes and flagged partial", {
  segs = seg_dt(c("chr1", "chr1", "chr1"), c(1, 1001, 2001), c(1000, 2000, 3000),
                c(1, 4, 1), c(1, 2, 0))
  p = list(lym = fake_panel(c("SPAN", "ONE", "GAP"), c("chr1", "chr1", "chr2"),
                            c(500, 100, 10), c(2500, 200, 20), c("SPAN", "ONE", "GAP"),
                            role_gene = "SPAN", role = "oncogene"))
  g = cna_gene_calls(segs, p, ploidy = 2, sex = "female", source = "wakhan")$genes
  setkey(g, gene)
  expect_equal(g["SPAN"]$call, "amplification")   # max_cn 6
  expect_equal(g["SPAN"]$min_cn, 1); expect_equal(g["SPAN"]$max_cn, 6)
  expect_equal(g["SPAN"]$n_segments, 3L)
  expect_true(g["SPAN"]$partial)
  expect_true(g["SPAN"]$loh)                       # the 1/0 segment
  expect_equal(g["SPAN"]$role, "oncogene")
  expect_equal(g["ONE"]$call, "neutral"); expect_false(g["ONE"]$partial)
  expect_true(is.na(g["ONE"]$role))
  # A gene in a gap keeps its row with no call
  expect_true(is.na(g["GAP"]$call)); expect_equal(g["GAP"]$n_segments, 0L)
  expect_equal(unique(g$source), "wakhan")
  expect_equal(names(g), names(.empty_cna_genes()))
})

test_that("ASCAT segments with bare contigs are chr-prefixed and uncapped", {
  s = parse_ascat_integer_segments(cna_fx_path("ascat", "S1.segments_raw.txt"))
  expect_equal(s$chrom, c("chr1", "chr1", "chr2", "chrX"))
  expect_equal(s$total, c(2, 6, 1, 1))          # not capped at 4 like parse_ascat_segments()
  expect_equal(s$minor, c(1, 2, 0, 0))
  expect_equal(names(s), names(.empty_cn_segments()))
  expect_null(parse_ascat_integer_segments(NULL))
  expect_null(parse_ascat_integer_segments(tempfile()))
  m = parse_ascat_metrics(cna_fx_path("ascat", "S1.metrics.txt"))
  expect_equal(m$wgd, 1); expect_equal(m$loh, 0.12); expect_equal(m$gi, 0.3)
})

test_that("integer_profile.bed takes its header from the LAST # line", {
  s = parse_wakhan_integer_profile(cna_fx_path("wakhan", "solution_rank_1", "integer_profile.bed"))
  expect_equal(nrow(s), 5)
  expect_equal(s$chrom[1:3], rep("chr1", 3))
  expect_equal(s$cn_a, c(1, 3, 0, 1, 1)); expect_equal(s$cn_b, c(1, 3, 0, 0, 0))
  expect_equal(s$total, c(2, 6, 0, 1, 1))
  expect_equal(s$conf[2], 0.81)
  expect_equal(s$sv_ids[2], "severus_DUP1;severus_DEL2")
  expect_true(is.na(s$sv_ids[3]))               # empty trailing column
  raw = read_hash_header_table(cna_fx_path("wakhan", "solution_rank_1", "integer_profile.bed"))
  expect_equal(names(raw)[1], "chr")
  expect_null(parse_wakhan_integer_profile(tempfile()))
})

test_that("Wakhan gene rows with implausible states are dropped and counted", {
  g = parse_wakhan_genes(cna_fx_path("wakhan", "solution_rank_1", "genes", "genes_copynumber_states.bed"),
                         CNA_MAX_PLAUSIBLE_CN)
  expect_equal(g$n_dropped, 1L)
  expect_equal(sort(g$genes$gene), c("OTHER", "PANELONLY"))
  expect_equal(g$genes[gene == "PANELONLY"]$total, 4)
})

test_that("WGD and LOH fraction are length-weighted over the autosomes", {
  # 60% of autosomal length with major >= 2, 25% LOH; chrX ignored
  segs = seg_dt(c("chr1", "chr2", "chr3", "chrX"), c(1, 1, 1, 1), c(600, 150, 250, 10000),
                c(2, 1, 1, 0), c(2, 0, 1, 0))
  g = cna_genome_metrics(segs)
  expect_true(g$wgd)
  expect_equal(g$major2_fraction, 0.6)
  expect_equal(g$loh_fraction, 0.15)
  segs2 = seg_dt(c("chr1", "chr2"), c(1, 1), c(400, 600), c(2, 1), c(1, 1))
  expect_false(cna_genome_metrics(segs2)$wgd)
  expect_true(is.na(cna_genome_metrics(NULL)$wgd))
})

test_that("Wakhan is preferred for segments and purity; every source is listed", {
  d = make_sample_dir("copy")
  inp = cna_locate(d, "S1")
  expect_match(inp$ascat_metrics, "S1\\.metrics\\.txt$")
  p = list(lym = fake_panel(c("FOCAL", "PANELONLY", "NOWHERE"), "chr1", 1200000, 1300000, "FOCAL"))
  x = cna_parse(inp, list(), cna_ctx(panels = p))
  expect_equal(x$source, "wakhan")
  expect_equal(x$purity, 0.57); expect_equal(x$ploidy, 2.02)
  expect_equal(x$purity_by_source$source, c("wakhan", "ascat", "savana"))
  expect_equal(x$purity_by_source$used, c(TRUE, FALSE, FALSE))
  expect_equal(x$purity_by_source[source == "savana"]$purity, 0.47)   # rank 1, not row 1
  expect_equal(x$purity_by_source[source == "ascat"]$purity, 0.81)
  expect_equal(nrow(x$segments), 5)
  expect_equal(x$genome$ascat_wgd, 1)
  expect_equal(x$wakhan_genes_dropped, 1L)
  setkey(x$genes, gene)
  expect_equal(x$genes["FOCAL"]$call, "amplification")
  # Symbol-only panel gene: called from Wakhan's gene file; a gene in neither has no row
  expect_equal(x$genes["PANELONLY"]$source, "wakhan-genes")
  expect_equal(x$genes["PANELONLY"]$call, "gain")
  expect_false("NOWHERE" %in% x$genes$gene)
  expect_length(x$profiles, 1)
  expect_equal(x$solutions$solution_rank, c(1L, 2L))
})

test_that("ASCAT supplies segments and purity when Wakhan is absent", {
  d = make_sample_dir("copy", wakhan = FALSE, savana = FALSE)
  inp = cna_locate(d, "S1")
  expect_true(is.na(inp$wakhan$dir))
  p = list(lym = fake_panel("FOCAL", "chr1", 1200000, 1300000, "FOCAL"))
  x = cna_parse(inp, list(), cna_ctx(panels = p))
  expect_equal(x$source, "ascat")
  expect_equal(x$purity, 0.81); expect_equal(x$ploidy, 2.1)
  expect_equal(x$purity_by_source$source, "ascat")
  expect_equal(x$genes$call, "amplification")
  expect_equal(x$genes$source, "ascat")
  expect_equal(nrow(x$wakhan_genes), 0)
  expect_length(x$profiles, 0)
})

test_that("nothing found returns NULL; panels without coordinates are reported", {
  d = cna_tempdir("empty", envir = environment())
  expect_null(cna_parse(cna_locate(d, "S1"), list(), cna_ctx()))
  d2 = make_sample_dir("copy", wakhan = FALSE, savana = FALSE)
  x = cna_parse(cna_locate(d2, "S1"), list(), cna_ctx(panels = list(sym = fake_panel(c("A", "B")))))
  expect_equal(x$panels_without_coords, "sym")
  expect_equal(nrow(x$genes), 0)
})

test_that("discovery works in the pipeline-staged layout of individually symlinked files", {
  d = make_sample_dir("staged")
  expect_true(nzchar(Sys.readlink(file.path(d, "wakhan", "solution_rank_1", "integer_profile.bed"))))
  inp = cna_locate(d, "S1")
  expect_match(inp$wakhan$profiles[[1]]$bed, "solution_rank_1/integer_profile.bed$")
  expect_match(inp$wakhan$profiles[[1]]$genes, "genes_copynumber_states.bed$")
  expect_true(is.na(inp$wakhan$heatmap))
  expect_match(inp$ascat_segments, "S1\\.segments_raw\\.txt$")
  expect_match(inp$ascat_metrics, "S1\\.metrics\\.txt$")
  x = cna_parse(inp, list(), cna_ctx())
  expect_equal(x$source, "wakhan")
  expect_equal(nrow(x$segments), 5)
  # The rank-2 solution is listed but has no profile in the staged layout
  expect_length(x$profiles, 1)
})

test_that("a symlinked solution_rank_1 directory is found by literal lookup", {
  d = make_sample_dir("copy")
  expect_true(nzchar(Sys.readlink(file.path(d, "wakhan", "solution_rank_1"))))
  hit = locate_wakhan_solution_file(file.path(d, "wakhan"), 1, "solution_2.02_0.57_0.92",
                                    "integer_profile.bed")
  expect_match(hit, "solution_rank_1/integer_profile.bed$")
  # Falls back to repository_name when the rank link is missing
  unlink(file.path(d, "wakhan", "solution_rank_1"))
  hit = locate_wakhan_solution_file(file.path(d, "wakhan"), 1, "solution_2.02_0.57_0.92",
                                    "integer_profile.bed")
  expect_match(hit, "solution_2.02_0.57_0.92/integer_profile.bed$")
  expect_null(locate_wakhan_solution_file(file.path(d, "wakhan"), 2, NA, "integer_profile.bed"))
})

test_that("Wakhan segments without a solutions table get an estimated ploidy", {
  d = make_sample_dir("staged", ascat = FALSE, savana = FALSE)
  unlink(file.path(d, "wakhan", "solutions_ranks.tsv"))
  x = cna_parse(cna_locate(d, "S1"), list(), cna_ctx())
  expect_equal(x$source, "wakhan")
  expect_true(x$ploidy_estimated)
  expect_true(is.na(x$purity))
  expect_equal(x$ploidy, round((2 * 1e6 + 6 * 1e6 + 0 + 1 * 5e6) / 8e6, 2))
})

test_that("the profile plot builds from segments and is NULL when there are none", {
  s = parse_wakhan_integer_profile(cna_fx_path("wakhan", "solution_rank_1", "integer_profile.bed"))
  p = cna_profile_plot(s, chromosomes_for_sex("female"),
                       setNames(rep(5e6, 23), c(paste0("chr", 1:22), "chrX")))
  expect_s3_class(p, "ggplot")
  expect_null(cna_profile_plot(.empty_cn_segments(), "chr1", c(chr1 = 1e6)))
  expect_null(cna_profile_plot(NULL, "chr1", c(chr1 = 1e6)))
})

test_that("the circos copy-number ring follows the CN section's source (REP-8)", {
  wk = list(source = "wakhan", segments = data.table(
    chrom = c("chr1", "chr2", "chr3"), start = c(1, 1, 1), end = c(100, 200, 300),
    cn_a = c(1, 5, NA), cn_b = c(3, 2, NA), total = c(4, 7, NA), minor = c(1, 2, NA)))
  r = circos_cn_segments(wk, ascat_segments = data.table(chr = "chr9"))
  expect_equal(r$source, "wakhan")
  expect_equal(r$segments$chr, c("chr1", "chr2"))           # a segment without a total is not drawn
  expect_equal(r$segments$major_cn, c(3, 4))                 # higher haplotype, capped at 4
  expect_equal(r$segments$minor_cn, c(1, 2))
  expect_equal(r$segments$total_cn, c(4, 4))
  # skip_ascat: Wakhan still draws a ring; no CN section at all: ASCAT's own segments, else nothing
  expect_equal(circos_cn_segments(NULL, data.table(chr = "chr9"))$source, "ascat")
  expect_null(circos_cn_segments(NULL, NULL))
  expect_null(circos_cn_segments(list(source = NA_character_, segments = data.table()), NULL))
})

test_that("a symbol-only list takes gene spans from another list before Wakhan's file (REP-9)", {
  segs = seg_dt(c("chr9", "chr3"), c(21e6, 160e6), c(23e6, 190e6), c(0, 2), c(0, 1))
  builtin = fake_panel(c("CDKN2A", "PIK3CA"), c("chr9", "chr3"), c(21967752, 179148114),
                       c(21995301, 179240093), c("CDKN2A", "PIK3CA"))
  custom = fake_panel(c("cdkn2a", "PIK3CA", "ONLYWAKHAN", "NOWHERE"))
  wg = data.table(chrom = "chr5", start = 1e6, end = 2e7, gene = "ONLYWAKHAN",
                  hp1_state = 2, hp2_state = 2, total = 4, minor = 2)
  g = cna_gene_calls(segs, list(sarcoma = builtin, mine = custom), ploidy = 2, sex = "female",
                     source = "wakhan", wakhan_genes = wg)$genes[panel == "mine"]
  setkey(g, gene)
  # The list's own spelling, the builtin's span, called on the segments like any coordinate row
  expect_equal(g["cdkn2a"]$call, "homozygous deletion")
  expect_equal(g["cdkn2a"]$start, 21967752)
  expect_equal(g["cdkn2a"]$span_from, "sarcoma")
  expect_equal(g["PIK3CA"]$source, "wakhan")
  # Only what no list places falls back to Wakhan's file, labelled as a segment
  expect_equal(g["ONLYWAKHAN"]$source, "wakhan-genes")
  expect_equal(g["ONLYWAKHAN"]$span_from, "Wakhan segment")
  expect_false("NOWHERE" %in% g$gene)
  # A list's own coordinate rows carry no span_from
  b = cna_gene_calls(segs, list(sarcoma = builtin, mine = custom), 2, "female", "wakhan")$genes[panel == "sarcoma"]
  expect_true(all(is.na(b$span_from)))
})

test_that("the ASCAT tab says 'skipped' without claiming a failure", {
  none = ascat_missing_note(list(), "wakhan")
  expect_match(none, "ASCAT was skipped \\(--skip_ascat\\) or produced nothing")
  expect_match(none, "comes from Wakhan")
  expect_false(grepl("failed", none))
  expect_match(ascat_missing_note(list(ascat_segments = "x.segments_raw.txt"), "ascat"),
               "^ASCAT ran .* plots were not found\\.$")
})

test_that("purity disagreement compares each caller's top solution only", {
  # Fixtures: Wakhan rank 1 0.57 (rank 2 is 0.99), ASCAT 0.81, SAVANA rank 1 0.47 (row 1 is 0.52)
  x = cna_parse(cna_locate(make_sample_dir("staged"), "S1"), list(), cna_ctx())
  o = build_overview(NULL, list(cna = x), list(), NULL, NULL, list(), list())
  pc = o$purity_check
  expect_equal(pc$sources$source, c("wakhan", "ascat", "savana"))
  expect_equal(pc$sources$purity, c(0.57, 0.81, 0.47))
  expect_equal(pc$spread, 0.34)
  expect_true(pc$flag)
  expect_equal(pc$used, "wakhan")
  expect_equal(o$biomarkers$purity$status, "warn")
  # SAVANA and ASCAT skipped: one caller, no check
  x = cna_parse(cna_locate(make_sample_dir("staged", ascat = FALSE, savana = FALSE), "S1"),
                list(), cna_ctx())
  expect_null(build_overview(NULL, list(cna = x), list(), NULL, NULL, list(), list())$purity_check)
})
