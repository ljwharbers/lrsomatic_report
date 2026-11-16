# WP6: slim small-variant embedding (slim_variant_table, clinvar_is_pathogenic,
# panel_snv_genes) and the circos density binning (bin_snv_density, circos_density_cap)

slim_fixture = function() {
  vt = data.table(
    symbol  = c("TP53", "KRAS", "FOO", "BAR", "MYC", "baz", "QUX", "ZED"),
    chrom   = "chr1",
    pos     = 1:8,
    ref     = "A", alt = "T",
    impact  = c("HIGH", "MODERATE", "LOW", "MODIFIER", "MODIFIER", "LOW", "LOW", NA),
    clinvar = c(NA, NA, "Pathogenic/Likely_pathogenic",
                "Conflicting_classifications_of_pathogenicity", NA, NA,
                "benign,likely_pathogenic", NA)
  )
  setattr(vt, "vep_fields", c(clinvar = "ClinVar_CLNSIG"))
  vt
}

# A panel as it comes back from the Quarto YAML round trip: plain list, vectors as lists
yaml_panel = function(genes, genes_snv = genes, has_scopes = FALSE) {
  list(name = "p", has_coords = FALSE, has_scopes = has_scopes,
       genes = as.list(genes), genes_snv = if (length(genes_snv)) as.list(genes_snv) else NULL)
}

test_that("clinvar_is_pathogenic matches tokens, not substrings", {
  x = c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
        "likely_pathogenic", "pathogenic,risk_factor", "Pathogenic,_low_penetrance",
        "Uncertain_significance,Pathogenic", "Pathogenic|other",
        "Conflicting_classifications_of_pathogenicity",
        "Conflicting_interpretations_of_pathogenicity", "Benign/Likely_benign",
        "Uncertain_significance", "", NA)
  expect_equal(clinvar_is_pathogenic(x),
               c(rep(TRUE, 8), rep(FALSE, 6)))
  expect_identical(clinvar_is_pathogenic(character(0)), logical(0))
  expect_identical(clinvar_is_pathogenic(NA), FALSE)
})

test_that("slim keeps exactly HIGH/MODERATE, ClinVar pathogenic and panel-gene rows", {
  vt = slim_fixture()
  out = slim_variant_table(vt, list(lymphoid = yaml_panel(c("myc", "BAZ"))), "slim")
  # TP53/KRAS by impact, FOO/QUX by ClinVar, MYC/baz by panel (case-insensitive, any impact)
  expect_equal(out$symbol, c("TP53", "KRAS", "FOO", "MYC", "baz", "QUX"))
  expect_false("BAR" %in% out$symbol)      # "Conflicting...pathogenicity" is not pathogenic
  expect_false("ZED" %in% out$symbol)      # NA impact, nothing else
  expect_equal(names(out), names(vt))      # columns never dropped or reordered
})

test_that("panel-gene rows are kept regardless of impact, and only SNV-scope genes count", {
  vt = slim_fixture()
  no_panel = slim_variant_table(vt, list(), "slim")
  expect_false(any(c("MYC", "baz") %in% no_panel$symbol))

  # Scoped panel: MYC is SV-only (absent from genes_snv), BAZ is SNV-scope
  scoped = yaml_panel(c("MYC", "BAZ"), genes_snv = "BAZ", has_scopes = TRUE)
  out = slim_variant_table(vt, list(p = scoped), "slim")
  expect_true("baz" %in% out$symbol)
  expect_false("MYC" %in% out$symbol)

  # Scoped panel whose SNV set is empty (NULL after YAML): no fallback to `genes`
  sv_only = yaml_panel(c("MYC", "BAZ"), genes_snv = character(0), has_scopes = TRUE)
  expect_false(any(c("MYC", "baz") %in% slim_variant_table(vt, list(p = sv_only))$symbol))

  # Unscoped panel lacking genes_snv entirely falls back to `genes`
  legacy = list(name = "old", genes = list("MYC"))
  expect_true("MYC" %in% slim_variant_table(vt, list(old = legacy))$symbol)
})

test_that("slim preserves attributes and records n_full / n_embedded", {
  vt = slim_fixture()
  setattr(vt, "some_other", "x")
  out = slim_variant_table(vt, list(), "slim")
  expect_equal(attr(out, "vep_fields"), c(clinvar = "ClinVar_CLNSIG"))
  expect_equal(attr(out, "some_other"), "x")
  expect_equal(attr(out, "n_full"), 8L)
  expect_equal(attr(out, "n_embedded"), nrow(out))
  expect_equal(nrow(out), 4L)
  expect_true(is.data.table(out))
})

test_that("full mode is the identity (plus attributes)", {
  vt = slim_fixture()
  out = slim_variant_table(vt, list(), "full")
  expect_identical(lapply(out, identity), lapply(vt, identity))
  expect_equal(attr(out, "vep_fields"), attr(vt, "vep_fields"))
  expect_equal(attr(out, "n_full"), 8L)
  expect_equal(attr(out, "n_embedded"), 8L)
  expect_error(slim_variant_table(vt, list(), "bogus"), "unknown mode")
})

test_that("the input is not modified by reference", {
  vt = slim_fixture()
  before = copy(vt)
  out = slim_variant_table(vt, list(p = yaml_panel("MYC")), "slim")
  out[, symbol := "CHANGED"]
  out_full = slim_variant_table(vt, list(), "full")
  out_full[, impact := "CHANGED"]
  expect_equal(vt, before)
  expect_null(attr(vt, "n_full"))
  expect_null(attr(vt, "n_embedded"))
})

test_that("NULL, empty and annotation-free inputs", {
  expect_null(slim_variant_table(NULL, list()))
  empty = slim_variant_table(data.table(), list())
  expect_equal(nrow(empty), 0L)
  expect_equal(attr(empty, "n_full"), 0L)
  expect_equal(attr(empty, "n_embedded"), 0L)
  # No clinvar column (plugin not run) and NULL panels: impact alone decides
  vt = slim_fixture()[, clinvar := NULL]
  expect_equal(slim_variant_table(vt, NULL)$symbol, c("TP53", "KRAS"))
})

test_that("write_full_variant_table round-trips the full table through fread", {
  vt = slim_fixture()
  vt[, vaf := c(0.1, 0.25, NA, 0.5, 1, 0, 0.333, NA)]
  f = tempfile(fileext = ".tsv.gz")
  on.exit(unlink(f))
  write_full_variant_table(vt, f)
  back = fread(f, na.strings = c("", "NA"))
  expect_equal(nrow(back), nrow(vt))
  expect_equal(names(back), names(vt))
  expect_equal(back$symbol, vt$symbol)
  expect_equal(back$clinvar, vt$clinvar)
  expect_equal(back$vaf, vt$vaf)
})

# ---- circos density binning -------------------------------------------------------

test_that("bin_snv_density counts SNVs per window, including empty windows", {
  lens = c(chr1 = 2500, chr2 = 1000)
  dt = data.table(chrom = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr3", "chr1"),
                  pos   = c(1, 1000, 1001, 2500, 1500, 999, 5, 10),
                  ref   = c("A", "C", "G", "T", "A", "C", "A", "AT"),
                  alt   = c("G", "T", "A", "C", "G", "T", "G", "A"))
  # duplicated record is counted once
  dt = rbind(dt, dt[1])
  b = bin_snv_density(dt, lens, bin = 1000)
  expect_equal(names(b), c("chrom", "start", "end", "n", "per_mb"))
  expect_equal(b$chrom, c("chr1", "chr1", "chr1", "chr2"))
  expect_equal(b$start, c(0, 1000, 2000, 0))
  expect_equal(b$end,   c(1000, 2000, 2500, 1000))
  # chr1: pos 1, 1000 -> win 0; 1001, 1500 -> win 1; 2500 -> win 2; indel and chr3 ignored
  expect_equal(b$n, c(2L, 2L, 1L, 1L))
  expect_equal(b$per_mb, b$n * 1e6 / 1000)
})

test_that("bin_snv_density handles empty input and follows chrom_lens order", {
  lens = c(chrX = 3000, chr1 = 1000)
  b = bin_snv_density(NULL, lens, bin = 1000)
  expect_equal(unique(b$chrom), c("chrX", "chr1"))
  expect_true(all(b$n == 0L))
  expect_equal(nrow(bin_snv_density(NULL, NULL)), 0L)
  # Default window is CIRCOS_SNV_BIN: ~3,100 windows for a whole genome, whatever the load
  expect_equal(nrow(bin_snv_density(NULL, c(chr1 = 248e6))), 248L)
})

test_that("circos_density_cap clips outliers and never drops below 1", {
  expect_equal(circos_density_cap(c(rep(5, 200), 300)), 5)
  expect_equal(circos_density_cap(numeric(0)), 1)
  expect_equal(circos_density_cap(c(0, 0, 0)), 1)
  expect_equal(circos_density_cap(c(0.2, 0.4)), 1)
  # Zeros (empty windows) do not drag the cap down
  expect_equal(circos_density_cap(c(rep(0, 1000), rep(8, 50))), 8)
})
