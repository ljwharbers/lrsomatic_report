# Content checks on the builtin gene panels in assets/gene_lists/. The loader's own
# behaviour is tested in test-utils.R; these pin what the shipped files contain.

builtin_panel_files = function() {
  f = Sys.glob(file.path(gene_lists_root, "*.tsv"))
  ref = vapply(f, function(x) .split_panel_filename(x)$reference, character(1))
  data.table(path = f, name = vapply(f, function(x) .split_panel_filename(x)$name, character(1)),
             reference = ref)
}

# The raw table, comment header skipped, every column as character
read_panel_raw = function(path) {
  fread(path, skip = .panel_header(path)$n_comment, sep = "\t", colClasses = "character",
        na.strings = NULL)
}

ref_chroms = function(ref) {
  cl = fread(file.path(repo_root, "assets", "references", ref, "chrom_lengths.tsv"),
             header = FALSE, col.names = c("chrom", "length"))
  setNames(cl$length, cl$chrom)
}

test_that("every builtin panel file is reference-suffixed and loads on its own reference", {
  files = builtin_panel_files()
  expect_gt(nrow(files), 0)
  expect_false(anyNA(files$reference))
  for (i in seq_len(nrow(files))) {
    # load_gene_panel() directly: load_all_gene_panels() turns a load error into a message
    expect_no_error(load_gene_panel(files$path[i], files$reference[i]))
  }
  for (ref in c("hg38", "t2t")) {
    panels = load_all_gene_panels(gene_lists_root, ref)
    expect_setequal(names(panels), unique(files$name))
    for (p in panels) expect_equal(p$reference, ref)
  }
})

test_that("builtin panels have no duplicated gene rows", {
  files = builtin_panel_files()
  for (f in files$path) {
    g = read_panel_raw(f)$gene
    expect_false(anyDuplicated(g) > 0, info = paste(basename(f), ":", g[duplicated(g)]))
  }
})

test_that("builtin panel roles are within the allowed set", {
  files = builtin_panel_files()
  for (i in seq_len(nrow(files))) {
    raw = read_panel_raw(files$path[i])
    if ("role" %in% names(raw))
      expect_true(all(raw$role %in% c("oncogene", "tsg", "both", "")), info = files$path[i])
    r = panel_roles(load_gene_panel(files$path[i], files$reference[i]))
    expect_true(all(r$role %in% c("oncogene", "tsg", "both")), info = files$path[i])
  }
})

test_that("the lymphoid role column is the twist class from notes, and notes are untouched", {
  for (ref in c("hg38", "t2t")) {
    raw = read_panel_raw(file.path(gene_lists_root, sprintf("lymphoid.%s.tsv", ref)))
    cls = sub(";.*$", "", raw$notes)
    want = unname(c(OG = "oncogene", TSG = "tsg", `TSG/OG` = "both")[cls])
    want[is.na(want)] = ""
    expect_equal(raw$role, want)
    expect_equal(as.vector(table(factor(raw$role, c("oncogene", "tsg", "both", "")))),
                 c(74L, 64L, 10L, 86L))
    p = load_gene_panel(file.path(gene_lists_root, sprintf("lymphoid.%s.tsv", ref)), ref)
    expect_true(isTRUE(p$has_roles))
    roles = panel_roles(p)
    expect_equal(roles[gene == "MYD88", role], "oncogene")
    expect_equal(roles[gene == "TP53",  role], "tsg")
    expect_equal(roles[gene == "RHOA",  role], "both")
    # SV-only rows carry no role: none was given in the source
    expect_false("IGH" %in% roles$gene)
  }
})

test_that("builtin panel coordinates are sane for their reference", {
  files = builtin_panel_files()
  for (i in seq_len(nrow(files))) {
    lens = ref_chroms(files$reference[i])
    iv = panel_gene_intervals(load_gene_panel(files$path[i], files$reference[i]), "all")
    info = basename(files$path[i])
    expect_true(all(iv$chrom %in% names(lens)), info = info)
    expect_true(all(iv$start >= 1L & iv$start < iv$end), info = info)
    expect_true(all(iv$end <= lens[iv$chrom]), info = info)
  }
})

test_that("lymphoid snv rows now carry coordinates on both references", {
  for (ref in c("hg38", "t2t")) {
    p = load_gene_panel(file.path(gene_lists_root, sprintf("lymphoid.%s.tsv", ref)), ref)
    all_iv = panel_gene_intervals(p, "all")
    # 128 before (SV-scoped rows only); every row has coordinates now
    expect_gt(nrow(all_iv), 128L)
    expect_equal(nrow(all_iv), 234L)
    expect_true(all(c("MYD88", "NOTCH1", "PRKCB", "U2AF1", "CRLF2") %in% all_iv$gene))
    # CRLF2 is pseudoautosomal; the chrX copy is the one shipped
    expect_equal(all_iv[gene == "CRLF2", chrom], "chrX")
  }
  # the same gene sits on the same chromosome in both references
  h = panel_gene_intervals(load_gene_panel(file.path(gene_lists_root, "lymphoid.hg38.tsv"), "hg38"), "all")
  t = panel_gene_intervals(load_gene_panel(file.path(gene_lists_root, "lymphoid.t2t.tsv"), "t2t"), "all")
  m = merge(h, t, by = "gene")
  expect_equal(nrow(m), 234L)
  expect_equal(m$chrom.x, m$chrom.y)
})

test_that("coordinates on snv rows do not change SV matching", {
  # SV intervals before this change: lymphoid 128 on each reference, sarcoma 140
  before = c(lymphoid = 128L, sarcoma = 140L)
  for (ref in c("hg38", "t2t")) {
    for (nm in names(before)) {
      p = load_gene_panel(file.path(gene_lists_root, sprintf("%s.%s.tsv", nm, ref)), ref)
      sv_iv = panel_gene_intervals(p, "sv")
      expect_equal(nrow(sv_iv), before[[nm]], info = paste(nm, ref))
      expect_equal(nrow(panel_intervals(p)), before[[nm]], info = paste(nm, ref))
      expect_setequal(sv_iv$gene, p$genes_sv)
    }
    p = load_gene_panel(file.path(gene_lists_root, sprintf("lymphoid.%s.tsv", ref)), ref)
    # MYD88 has a span now, but stays snv-scoped: placed on the genome, never SV-matched
    expect_true("MYD88" %in% panel_gene_intervals(p, "all")$gene)
    expect_false("MYD88" %in% panel_gene_intervals(p, "sv")$gene)
    expect_false("MYD88" %in% p$genes_sv)
  }

  # and in isolation: a coordinate-carrying snv row is not an SV interval
  tsv = tempfile(fileext = ".tsv")
  on.exit(unlink(tsv))
  writeLines(c("gene\tchrom\tstart\tend\tapplies_to",
               "SVG\tchr1\t100\t200\tsv",
               "SNVG\tchr2\t300\t400\tsnv"), tsv)
  p = load_gene_panel(tsv)
  expect_equal(panel_intervals(p)$gene, "SVG")
  expect_setequal(panel_gene_intervals(p, "all")$gene, c("SVG", "SNVG"))
})

test_that("the pipeline-repo corrections are present", {
  for (ref in c("hg38", "t2t")) {
    g = read_panel_raw(file.path(gene_lists_root, sprintf("lymphoid.%s.tsv", ref)))$gene
    expect_true(all(c("PRKCB", "DDX6", "TENT5C") %in% g))
    expect_false(any(c("PRKBC", "RCK", "FAM46C") %in% g))
  }
  # EWSR1 / KDM6B / SIK3 carried CHM13 coordinates in the hg38 file
  iv = panel_gene_intervals(load_gene_panel(file.path(gene_lists_root, "lymphoid.hg38.tsv"), "hg38"), "sv")
  expect_equal(iv[gene == "EWSR1", start], 29268009L)
  expect_equal(iv[gene == "KDM6B", start], 7834217L)
  expect_equal(iv[gene == "SIK3",  start], 116843402L)
})
