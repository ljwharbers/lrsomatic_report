# Genome instability: ShatterSeek (R/sections/shatterseek.R) and ecDNA (R/sections/ecdna.R).
# Fixtures in tests/fixtures/instability/ are hand-made to the formats the pipeline writes:
# assets/shatterseek/run_shatterseek.R and the CoRAL / AmpliconClassifier 2.0.0 modules.

instab_fx = file.path(repo_root, "tests", "fixtures", "instability")

run_section = function(id, sample_dir, sample_id = "S1") {
  s = SECTIONS[[id]]
  run_section_parse(s, s$locate(sample_dir, sample_id), list(), list(embed = "slim"))
}

# Mirror a fixture tree with every file an individual symlink (how the pipeline stages it)
symlink_tree = function(src) {
  dst = tempfile("staged")
  files = list.files(src, recursive = TRUE, all.files = TRUE)
  for (f in files) {
    dir.create(dirname(file.path(dst, f)), recursive = TRUE, showWarnings = FALSE)
    file.symlink(normalizePath(file.path(src, f)), file.path(dst, f))
  }
  dst
}

# --- ShatterSeek -----------------------------------------------------------------------

test_that("absent ShatterSeek output is not_run", {
  x = run_section("shatterseek", file.path(instab_fx, "ec_detected"))
  expect_identical(x$status, "not_run")
  expect_match(shatterseek_verdict(x), "not run")
})

test_that("one 23-row ShatterSeek file parses to canonical columns with two high calls", {
  p = parse_shatterseek_file(file.path(instab_fx, "ss_two_combos/shatterseek/severus_ascat/S1.shatterseek.tsv"))
  expect_identical(p$combo, "severus_ascat")
  expect_equal(nrow(p$dt), 23)
  expect_true(all(c("chrom", "start", "end", "n_sv", "n_tra", "n_oscillating", "hc1", "hc2", "hc3",
                    "criteria", "confidence", "q_fragment_joins") %in% names(p$dt)))
  expect_identical(p$dt$chrom[1:3], c("chr1", "chr2", "chr3"))
  expect_identical(p$dt[confidence == "high", chrom], c("chr3", "chr8"))
  expect_identical(p$dt[chrom == "chr8", criteria], "HC1,HC2")
  expect_true(p$dt[chrom == "chr8", hc2])
  expect_true(is.na(p$dt[chrom == "chr21", start]))   # chromosome without a cluster
})

test_that("two combos that disagree are both kept, labelled, and flagged", {
  x = run_section("shatterseek", file.path(instab_fx, "ss_two_combos"))
  expect_identical(x$status, "ok")
  expect_identical(names(x$combos), c("severus_ascat", "severus_wakhan"))
  expect_identical(unname(x$labels), c("Severus + ASCAT", "Severus + Wakhan"))
  expect_true(all(c("combo", "chrom", "confidence") %in% names(x$high)))
  expect_equal(nrow(x$high), 3)
  expect_identical(x$high[combo == "severus_wakhan", chrom], "chr3")
  expect_true(shatterseek_disagree(x))
  v = shatterseek_verdict(x)
  # The verdict names the PCAWG tier(s) each chromosome meets, not only "high"
  expect_match(v, "^High-confidence chromothripsis on chr3 \\[HC1\\] \\(all combinations\\)")
  expect_match(v, "chr8 \\[HC1,HC2\\] \\(Severus \\+ ASCAT only\\)")
})

test_that("renamed, missing and extra columns are tolerated; directory names the combo", {
  x = run_section("shatterseek", file.path(instab_fx, "ss_renamed"))
  expect_identical(x$status, "ok")
  expect_identical(names(x$combos), "savana")
  d = x$combos$savana
  expect_equal(nrow(d), 23)
  expect_false("my_extra_column" %in% names(d))
  expect_true(all(is.na(d$q_exp_cluster)))            # missing optional column -> NA
  expect_equal(d[chrom == "chr5", n_sv], 14)           # cluster_size spelling
  expect_equal(d[chrom == "chr5", n_oscillating], 9)   # case-insensitive match
  # no confidence column: derived from the tiers / criteria
  expect_identical(x$high$chrom, "chr5")
  expect_identical(x$high$confidence, "high")
  expect_identical(x$high$criteria, "HC3")
  expect_false(shatterseek_disagree(x))
  expect_identical(shatterseek_verdict(x), "High-confidence chromothripsis on chr5 [HC3] (SAVANA).")
})

test_that("a file without a chromosome column is skipped with a message", {
  expect_message(x <- run_section("shatterseek", file.path(instab_fx, "ss_nochrom")),
                 "no chromosome column")
  expect_identical(x$status, "not_run")
})

test_that("non-'none' tiers other than 'high' count as calls", {
  f = tempfile(fileext = ".shatterseek.tsv")
  writeLines(c("sample\tcn_source\tsv_source\tchrom\tchromothripsis_confidence",
               "S1\tascat\tseverus\tchr1\tlow", "S1\tascat\tseverus\tchr2\tnone",
               "S1\tascat\tseverus\tchr3\t"), f)
  p = parse_shatterseek_file(f)
  h = shatterseek_calls(list(severus_ascat = p$dt))
  expect_identical(h$chrom, "chr1")
  expect_identical(h$confidence, "low")
  x = list(status = "ok", combos = list(severus_ascat = p$dt), labels = c(severus_ascat = "Severus + ASCAT"),
           high = h)
  expect_match(shatterseek_verdict(x), "^Chromothripsis called on chr1 \\[low\\]")
})

test_that("no calls gives a 'No chromothripsis' verdict", {
  f = file.path(tempfile(), "shatterseek", "severus_ascat", "S1.shatterseek.tsv")
  dir.create(dirname(f), recursive = TRUE)
  d = fread(file.path(instab_fx, "ss_two_combos/shatterseek/severus_ascat/S1.shatterseek.tsv"))
  d[, chromothripsis_confidence := "none"]
  fwrite(d, f, sep = "\t")
  x = run_section("shatterseek", dirname(dirname(dirname(f))))
  expect_equal(nrow(x$high), 0)
  n = nrow(d)
  expect_identical(shatterseek_verdict(x),
                   sprintf("No chromothripsis called (Severus + ASCAT; %d chromosomes tested).", n))

  # Several combinations, all negative: one sentence, not one per combination
  root = tempfile()
  neg = list()
  for (cb in c("severus_ascat", "severus_wakhan")) {
    g = file.path(root, "shatterseek", cb, "S1.shatterseek.tsv")
    dir.create(dirname(g), recursive = TRUE)
    neg[[cb]] = fread(file.path(instab_fx, "ss_two_combos/shatterseek", cb, "S1.shatterseek.tsv"))
    neg[[cb]][, chromothripsis_confidence := "none"]
    fwrite(neg[[cb]], g, sep = "\t")
  }
  x = run_section("shatterseek", root)
  expect_identical(shatterseek_verdict(x),
                   sprintf("No chromothripsis called by any of the 2 caller combinations (%s; %d chromosomes tested each).",
                           paste(unname(x$labels), collapse = ", "), n))
  # Different numbers tested: one count per combination
  fwrite(neg$severus_wakhan[-1], file.path(root, "shatterseek", "severus_wakhan", "S1.shatterseek.tsv"),
         sep = "\t")
  x = run_section("shatterseek", root)
  expect_match(shatterseek_verdict(x), sprintf("%s: %d chromosomes tested", x$labels[["severus_wakhan"]], n - 1L),
               fixed = TRUE)
})

test_that("ShatterSeek is found through individually symlinked files", {
  x = run_section("shatterseek", symlink_tree(file.path(instab_fx, "ss_two_combos")))
  expect_identical(x$status, "ok")
  expect_length(x$combos, 2)
})

# --- ecDNA -----------------------------------------------------------------------------

test_that("no ecdna output is not_run", {
  x = run_section("ecdna", file.path(instab_fx, "ss_two_combos"))
  expect_identical(x$status, "not_run")
  expect_equal(nrow(x$amplicons), 0)
})

test_that("a 0-byte seed BED is none_detected, not not_run", {
  x = run_section("ecdna", file.path(instab_fx, "ec_empty_seeds"))
  expect_identical(x$status, "none_detected")
  expect_true(x$seeds_empty)
  expect_match(ecdna_verdict(x), "seed threshold")
})

test_that("classifier output builds the amplicon table; No amp/Invalid is dropped", {
  x = run_section("ecdna", file.path(instab_fx, "ec_detected"))
  expect_identical(x$status, "detected")
  a = x$amplicons
  expect_true(all(c("amplicon", "class", "ecdna", "bfb", "location", "oncogenes", "max_cn") %in% names(a)))
  expect_identical(a$amplicon, 1:2)
  expect_identical(a$class, c("ecDNA", "Linear"))
  expect_identical(a$ecdna, c(TRUE, FALSE))
  expect_identical(a$bfb, c(FALSE, FALSE))
  expect_identical(a$oncogenes, c("MYC", "MDM2"))
  expect_identical(a$location[1], "chr8:126500000-128200000, chr8:129000000-129400000")
  expect_equal(a$max_cn, c(32.4, 6.3))
  expect_equal(x$n_ecdna, 1)
  expect_equal(x$n_bfb, 0)
  expect_equal(x$n_invalid, 1)
  expect_identical(x$genes$gene, c("MYC", "PVT1", "MDM2"))
  expect_identical(x$genes$oncogene, c(TRUE, FALSE, TRUE))
  expect_identical(basename(x$plots), c("S1_amplicon1_graph.png", "S1_amplicon1_cycles.png"))
  expect_identical(ecdna_verdict(x), "1 ecDNA+ of 2 amplicons.")
})

test_that("a header-only classifier file is none_detected", {
  x = run_section("ecdna", file.path(instab_fx, "ec_header_only"))
  expect_identical(x$status, "none_detected")
  expect_false(x$seeds_empty)
  expect_match(ecdna_verdict(x), "No focal amplicons classified")
})

test_that("amplicons with no ecDNA+/BFB+ say so; renamed columns resolve; CoRAL fills gaps", {
  x = run_section("ecdna", file.path(instab_fx, "ec_no_ecdna"))
  expect_identical(x$status, "detected")
  expect_equal(nrow(x$amplicons), 2)
  expect_equal(x$n_ecdna, 0)
  expect_identical(x$amplicons$class, c("Complex-non-cyclic", "Linear"))
  expect_identical(x$amplicons$location[2], "chr12:68800000-69300000")   # from CoRAL summary
  expect_equal(x$amplicons$max_cn, c(31.5, 6))
  expect_identical(ecdna_verdict(x), "2 amplicons found, none classified as ecDNA or BFB.")
})

test_that("CoRAL without a classifier is unclassified, with CoRAL's amplicons", {
  x = run_section("ecdna", file.path(instab_fx, "ec_coral_only"))
  expect_identical(x$status, "unclassified")
  expect_equal(nrow(x$amplicons), 2)
  expect_true(all(is.na(x$amplicons$ecdna)))
  expect_identical(x$amplicons$source, c("CoRAL", "CoRAL"))
  expect_equal(x$amplicons$size_bp, c(2100002, 500001))
  expect_true(is.na(x$n_ecdna))
  expect_match(ecdna_verdict(x), "not classified")
})

test_that("ecDNA is found through individually symlinked files", {
  x = run_section("ecdna", symlink_tree(file.path(instab_fx, "ec_detected")))
  expect_identical(x$status, "detected")
  expect_length(x$plots, 2)
  x = run_section("ecdna", symlink_tree(file.path(instab_fx, "ec_empty_seeds")))
  expect_identical(x$status, "none_detected")
})

test_that("generic _gene_list.tsv outside ecdna/ is not picked up", {
  d = tempfile("ecgen")
  dir.create(file.path(d, "other"), recursive = TRUE)
  writeLines("gene\nTP53", file.path(d, "other", "S1_gene_list.tsv"))
  expect_identical(run_section("ecdna", d)$status, "not_run")
})

test_that("real AITL6 (0-byte seeds) parses to none_detected", {
  aitl6 = "/gpfs1/project/prj_00009/home/averham/LR_SOMATIC_T2T/AITL6"
  skip_if_not(dir.exists(file.path(aitl6, "ecdna")), "AITL6 sample dir not available")
  s = SECTIONS$ecdna
  inputs = s$locate(file.path(aitl6, "ecdna"), "AITL6")
  x = s$parse(inputs, list(), list())
  expect_identical(x$status, "none_detected")
})

test_that("the best tier of a call and the minimum-tier views", {
  expect_identical(shatterseek_best_tier(c("HC1,HC2", "HC3", "HC2,HC3", "", NA)),
                   c(1L, 3L, 2L, NA, NA))
  expect_identical(shatterseek_tiers(c("HC3", "HC1,HC2")), "HC1,HC2,HC3")
  x = run_section("shatterseek", file.path(instab_fx, "ss_two_combos"))
  # chr3 (HC1) in both combinations, chr8 HC1,HC2 in one
  expect_identical(shatterseek_at_tier(x, 3L), x)
  expect_equal(nrow(shatterseek_at_tier(x, 1L)$high), 3)
  y = x
  y$high = copy(x$high)[chrom == "chr3", criteria := "HC3"]
  h1 = shatterseek_at_tier(y, 1L)
  expect_identical(h1$high$chrom, "chr8")
  expect_equal(h1$n_below_tier, 2L)
  expect_match(shatterseek_verdict(h1), "chr8 \\[HC1,HC2\\]")
  # No call left at this tier: say so, never "No chromothripsis called"
  y$high = y$high[chrom == "chr3"]
  v = shatterseek_verdict(shatterseek_at_tier(y, 2L))
  expect_match(v, "^No chromothripsis call meets HC1–HC2; 2 calls meet only a lower tier")
})
