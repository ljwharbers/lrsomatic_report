# Mutational signatures (R/sections/signatures.R); fixtures trimmed from pilot15 (P15)

sig_fixture = file.path(repo_root, "tests", "fixtures", "signatures", "p15")
sig_assets  = file.path(repo_root, "assets", "signatures")
sig_fit_dir = function(root, ctx) file.path(root, "signatures", "assignment", "COSMIC_v3.6", ctx,
                                            "Assignment_Solution")

# Copy the fixture sample dir so a test can delete or rewrite parts of it
sig_copy_fixture = function() {
  d = tempfile("sigfix")
  dir.create(d)
  file.copy(file.path(sig_fixture, "signatures"), d, recursive = TRUE)
  d
}

test_that("activities melt to a long table of non-zero signatures with pct", {
  act = parse_spa_activities(file.path(sig_fit_dir(sig_fixture, "SBS96"), "Activities",
                                       "Assignment_Solution_Activities.txt"))
  expect_equal(act$signature, c("SBS5", "SBS1"))
  expect_equal(act$mutations, c(20517, 2831))
  expect_equal(sum(act$pct), 100)
  expect_equal(act$pct[2], 100 * 2831 / 23348)
})

test_that("stats parse strips the % from the percent columns", {
  st = parse_spa_stats(file.path(sig_fit_dir(sig_fixture, "SBS96"), "Solution_Stats",
                                 "Assignment_Solution_Samples_Stats.txt"))
  expect_equal(st$total, 23348)
  expect_equal(st$cosine, 0.993)
  expect_equal(st$l1_pct, 10.944)
  expect_equal(st$l2_pct, 11.85)
  expect_equal(st$kl, 0.01035)
  expect_null(parse_spa_stats(file.path(tempdir(), "no_such_stats.txt")))
})

test_that("channel grouping yields 6/10/16 classes over 96/78/83 channels in COSMIC order", {
  for (k in names(SIG_CONTEXTS)) {
    ch = sig_channels(k)
    g  = sig_channel_group(ch, k)
    expect_false(anyNA(g))
    expect_equal(length(unique(g)), c(SBS = 6, DBS = 10, ID = 16)[[k]])
    # groups appear as contiguous runs, in the plotting order of the colour palette
    expect_equal(rle(g)$values, sig_group_levels(k))
  }
  expect_length(sig_channels("SBS96"), 96)
  expect_length(sig_channels("DBS78"), 78)
  expect_length(sig_channels("ID83"), 83)
  # substitution class first, not MatrixGenerator's alphabetical 5'-base-first order
  expect_equal(sig_channels("SBS")[1:3], c("A[C>A]A", "A[C>A]C", "A[C>A]G"))
  expect_equal(sig_channels("SBS")[17], "A[C>G]A")
  expect_equal(sig_channel_group(c("1:Del:C:0", "5:Ins:R:5", "5:Del:M:5", "3:Del:R:2"), "ID"),
               c("1bp Del C", "5+bp Ins repeat", "5+bp Del MH", "3bp Del repeat"))
  expect_equal(sig_channel_group(c("AC>CA", "TT>GG"), "DBS"), c("AC>NN", "TT>NN"))
  expect_true(is.na(sig_channel_group("garbage", "SBS")))
})

test_that("matrices parse into standard order and keep every count", {
  m = parse_sig_matrix(file.path(sig_fixture, "signatures/matrices/output/SBS/pilot15.SBS96.all"), "SBS")
  expect_equal(names(m), c("channel", "group", "count"))
  expect_equal(m$channel, sig_channels("SBS"))
  expect_equal(sum(m$count), 23348)
  expect_equal(m[channel == "A[C>A]A", count], 213)
  d = parse_sig_matrix(file.path(sig_fixture, "signatures/matrices/output/DBS/pilot15.DBS78.all"), "DBS78")
  expect_equal(nrow(d), 78)
  i = parse_sig_matrix(file.path(sig_fixture, "signatures/matrices/output/ID/pilot15.ID83.all"), "ID83")
  expect_equal(nrow(i), 83)
  expect_equal(i[channel == "1:Del:C:0", count], 77)
})

test_that("the full P15 fixture parses with the expected flags and dominant set", {
  inp = locate_signatures(sig_fixture, "pilot15")
  expect_equal(inp$version, "3.6")
  expect_setequal(names(inp$contexts), c("SBS", "DBS", "ID"))
  r = parse_signatures(inp, sig_assets)
  expect_equal(r$cosmic_version, "3.6")
  expect_equal(r$lookup_version, "3.6")
  expect_length(r$contexts$SBS$flags, 0)
  # DBS: 43 mutations, cosine 0.732 -> both flags; ID: cosine 0.893 -> poor fit only
  expect_setequal(r$contexts$DBS$flags, c("few_mutations", "poor_fit"))
  expect_equal(r$contexts$ID$flags, "poor_fit")
  # only the unflagged SBS context reaches the overview
  expect_equal(r$dominant$signature, c("SBS5", "SBS1"))
  expect_equal(unique(r$dominant$context), "SBS")
  expect_equal(r$dominant$aetiology[2], "Clock-like (5mC deamination)")
  expect_true(all(c("context", "signature", "pct", "aetiology", "artefact") %in% names(r$dominant)))
  expect_s3_class(plot_sig_profile(r$contexts$ID$matrix, "ID"), "ggplot")
})

test_that("a missing DBS directory drops that context only", {
  d = sig_copy_fixture()
  unlink(file.path(d, "signatures/assignment/COSMIC_v3.6/DBS78"), recursive = TRUE)
  r = parse_signatures(locate_signatures(d, "pilot15"), sig_assets)
  expect_setequal(names(r$contexts), c("SBS", "ID"))
})

test_that("no signature output at all gives NULL", {
  d = tempfile("sigempty"); dir.create(d)
  expect_null(SECTIONS$signatures$parse(SECTIONS$signatures$locate(d, "x"), list(),
                                        list(repo_dir = repo_root)))
})

test_that("a missing stats file degrades to no cosine, total from the matrix", {
  d = sig_copy_fixture()
  unlink(file.path(sig_fit_dir(d, "SBS96"), "Solution_Stats"), recursive = TRUE)
  r = parse_signatures(locate_signatures(d, "pilot15"), sig_assets)
  sbs = r$contexts$SBS
  expect_null(sbs$files$stats)
  expect_true(is.na(sbs$stats$cosine))
  expect_equal(sbs$stats$total, 23348)
  expect_length(sbs$flags, 0)
})

test_that("all-zero activities mean no signature assigned, not an error", {
  f = tempfile(fileext = ".txt")
  writeLines(c("Samples\tSBS1\tSBS5", "s1\t0\t0"), f)
  act = parse_spa_activities(f)
  expect_equal(nrow(act), 0)
  expect_equal(names(act), c("signature", "mutations", "pct"))
  ann = annotate_sig_activities(act, load_sig_aetiology(sig_assets, "3.6"))
  expect_equal(nrow(ann), 0)
  expect_equal(nrow(sig_dominant(list(SBS = list(activities = ann, flags = character())))), 0)
})

test_that("an unknown signature reads 'not in lookup'; artefacts are flagged", {
  lk  = load_sig_aetiology(sig_assets, "3.6")
  act = data.table(signature = c("SBS999z", "SBS27", "SBS1"), mutations = c(50, 30, 20),
                   pct = c(50, 30, 20))
  ann = annotate_sig_activities(act, lk)
  expect_equal(ann$aetiology[1], "not in lookup")
  expect_false(ann$in_lookup[1])
  expect_false(ann$artefact[1])
  expect_true(ann$artefact[2])
  expect_false(ann$artefact[3])
  # no lookup at all still annotates
  expect_equal(annotate_sig_activities(act, NULL)$aetiology, rep("not in lookup", 3))
  # an unbundled version falls back to the highest bundled lookup
  expect_equal(attr(suppressMessages(load_sig_aetiology(sig_assets, "3.4")), "version"), "3.6")
})

test_that("low-count and low-cosine thresholds raise their flags", {
  expect_equal(sig_context_flags("SBS", 199, 0.95), "few_mutations")
  expect_length(sig_context_flags("SBS", 200, 0.9), 0)
  expect_equal(sig_context_flags("ID", 5000, 0.85), "poor_fit")
  expect_setequal(sig_context_flags("DBS", 10, 0.5), c("few_mutations", "poor_fit"))
  expect_length(sig_context_flags("DBS", 10 * SIG_MIN_MUTATIONS[["DBS"]], NA), 0)
})

test_that("flagged contexts and sub-threshold signatures stay out of dominant", {
  mk = function(sig, pct, flags = character())
    list(activities = data.table(signature = sig, mutations = pct, pct = pct,
                                 aetiology = "x", artefact = FALSE), flags = flags)
  dom = sig_dominant(list(SBS = mk(c("SBS1", "SBS2", "SBS3"), c(90, 6, 4)),
                          ID  = mk("ID1", 100, "poor_fit"),
                          DBS = mk("DBS1", 100)))
  expect_equal(dom$context, c("SBS", "SBS", "DBS"))
  expect_equal(dom$signature, c("SBS1", "SBS2", "DBS1"))
})

test_that("the lookup covers every signature column in the fixture headers", {
  lk = load_sig_aetiology(sig_assets, "3.6")
  expect_false(anyDuplicated(lk$signature) > 0)
  expect_true(all(lk$context %in% names(SIG_CONTEXTS)))
  for (ctx in SIG_CONTEXTS) {
    f = file.path(sig_fit_dir(sig_fixture, ctx), "Activities", "Assignment_Solution_Activities.txt")
    sigs = names(fread(f, nrows = 0))[-1]
    expect_true(length(sigs) > 0)
    missing = setdiff(sigs, lk$signature)
    expect_equal(missing, character(0), info = ctx)
  }
})
