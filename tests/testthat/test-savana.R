# SAVANA as the second SV caller. The fixture is 20 real records cut from a P15 SAVANA run
# (INFO trimmed to the keys the parser reads): mate pairs of every BP_NOTATION on one
# chromosome, two interchromosomal pairs, three INS and one orphaned mate.
fixture_dir = file.path(dirname(dirname(getwd())), "tests", "fixtures", "sv")
savana_vcf  = file.path(fixture_dir, "sample.classified.somatic.vcf")
savana_vep  = file.path(fixture_dir, "sample_VEP_SAVANA.vcf.gz")

savana_records = function() parse_sv_vcf(savana_vcf, "savana")

test_that("MATEID collapses SAVANA's breakend pairs to one row each", {
  r = savana_records()
  expect_equal(attr(r, "n_records"), 20L)
  # 8 pairs + 3 INS + 1 orphan
  expect_equal(nrow(r), 12L)
  expect_equal(anyDuplicated(r$id), 0L)
  expect_equal(r[id == "ID_25860_1", id_b], "ID_25860_2")
  # INS records are named ID_n_1 too, but have no mate and are not stem-paired
  expect_true(all(is.na(r[svtype == "INS", id_b])))
  expect_equal(sum(r$svtype == "INS"), 3L)
})

test_that("a non-PASS SAVANA record is left out, as for Severus", {
  # SAVANA writes its non-PASS calls (PREDICTED_NOISE, FAIL) to classified.vcf only; the
  # somatic VCFs of the campaign runs are all PASS, so this only guards the rule
  v = readLines(savana_vcf)
  i = grep("\tID_26414_1\t", v)
  v[i] = sub("\tPASS\t", "\tPREDICTED_NOISE\t", v[i])
  f = tempfile(fileext = ".vcf"); writeLines(v, f)
  on.exit(unlink(f))
  r = parse_sv_vcf(f, "savana")
  expect_equal(nrow(r), 11L)
  expect_false("ID_26414_1" %in% r$id)
  expect_equal(attr(r, "n_failed"), 1L)
  expect_equal(attr(savana_records(), "n_failed"), 0L)
})

test_that("same-chromosome pairs are typed from BP_NOTATION and become spans", {
  r = savana_records()
  expect_equal(r[id == "ID_25860_1", svtype], "DEL")   # +-
  expect_equal(r[id == "ID_29721_1", svtype], "DUP")   # -+
  expect_equal(r[id == "ID_26916_1", svtype], "INV")   # ++
  expect_equal(r[id == "ID_26917_1", svtype], "INV")   # --
  typed = r[svtype %in% c("DEL", "DUP", "INV")]
  expect_true(all(typed$svtype_raw == "BND"))
  # The span-style svclass the equivalent Severus types get, not a junction class
  expect_equal(typed$svclass, typed$svtype)
  expect_false(any(typed$svclass %in% SV_JUNCTION_CLASSES))
  # Left to right, sized
  expect_true(all(typed$pos_b > typed$pos_a))
  expect_equal(r[id == "ID_25860_1", sv_len], 421L)
  expect_equal(r[id == "ID_25860_1", c(pos_a, pos_b)], c(635796L, 636217L))
})

test_that("interchromosomal pairs stay translocations and orphans stay untyped", {
  r = savana_records()
  tr = r[svclass == "translocation"]
  expect_setequal(tr$id, c("ID_27202_1", "ID_21930_1"))
  expect_true(all(tr$svtype == "BND"))
  # SVLEN=0 on a translocation is "no length", not zero
  expect_true(all(is.na(tr$sv_len)))
  expect_equal(r[id == "ID_27202_1", chrom_b], "chr22")
  # The orphan's mate was cut: it keeps its ALT partner locus, as a Severus orphan does,
  # and is not orientation-typed (only a pair is)
  o = r[id == "ID_29727_1"]
  expect_equal(o$svtype, "BND")
  expect_true(is.na(o$id_b))
  expect_equal(o$svclass, "intra-chr breakend")
  expect_equal(o$pos_b, 196213044L)
})

test_that("VAF is the first TUMOUR_AF value and support is TUMOUR_READ_SUPPORT", {
  r = savana_records()
  # TUMOUR_AF=0.538,1 on the kept _1 record: this breakpoint first
  expect_equal(r[id == "ID_26916_1", vaf], 0.538)
  expect_equal(r[id == "ID_26916_1", support], 21L)
  expect_equal(r[id == "ID_26414_1", vaf], 0.5)
  expect_equal(r[id == "ID_26414_1", sv_len], 32L)
})

test_that("typed pairs flow into the contiguous circos ring and out of the breakend circos", {
  r = savana_records()
  tr = severus_circos_tracks(r)
  expect_equal(sort(unique(tr$nontrans$svtype)), c("DEL", "DUP", "INS", "INV"))
  expect_false(anyNA(tr$nontrans$circos_col))
  t = build_sv_table_from_vep(r, NULL)
  links = bnd_links(t, paste0("chr", c(1:22, "X", "Y")))
  expect_setequal(links$id, c("ID_27202_1", "ID_21930_1", "ID_29727_1"))
})

test_that("the SAVANA VEP annotation joins on record ID", {
  t = build_sv_table_from_vep(savana_records(), savana_vep)
  expect_equal(names(t), SV_TABLE_COLS)
  expect_equal(t[id == "ID_26555_1", gene_a], "KAZN")
  expect_equal(t[id == "ID_24263_1", impact], "MODIFIER")
  # ID_26588 is absent from the VEP file (as in the real one): kept, unannotated
  expect_true(is.na(t[id == "ID_26588_1", gene_a]))
  expect_equal(nrow(t), 12L)
  expect_gt(attr(t, "n_annotated"), 0L)
})

test_that("typed spans use the span window in panel matching, not the breakend window", {
  t = build_sv_table_from_vep(savana_records(), NULL)
  # A gene 500 kb from the DEL at chr1:635,796-636,217: inside a 1 Mb BND window, outside 100 kb
  p = list(has_coords = TRUE, genes = "FAR", interval_gene = "FAR",
           chrom = "chr1", start = 1136217L, end = 1137000L)
  expect_false(any(nzchar(sv_panel_hits(t[id == "ID_25860_1"], p))))
})

# ---- Discovery ------------------------------------------------------------

stage_sample = function(files) {
  d = tempfile("sample_")
  for (rel in names(files)) {
    dir.create(dirname(file.path(d, rel)), recursive = TRUE, showWarnings = FALSE)
    file.copy(files[[rel]], file.path(d, rel))
  }
  d
}

test_that("locate prefers the PoN-filtered SAVANA VCF and keeps the VEP files apart", {
  d = stage_sample(c(
    "variants/savana/sample.classified.somatic.vcf"              = savana_vcf,
    "variants/savana/sample.pon_filtered.classified.somatic.vcf" = savana_vcf,
    "variants/savana/sample.classified.vcf"                      = savana_vcf,
    "variants/savana/sample.classified.germline.vcf"             = savana_vcf,
    "variants/savana/sample.sv_breakpoints.vcf"                  = savana_vcf,
    "vep/SVs/sample_VEP_SAVANA.vcf.gz"                           = savana_vep
  ))
  on.exit(unlink(d, recursive = TRUE))
  inp = SECTIONS$sv$locate(d, "sample")
  expect_match(inp$callers$savana$vcf, "sample\\.pon_filtered\\.classified\\.somatic\\.vcf$")
  expect_true(inp$callers$savana$pon_filtered)
  expect_match(inp$callers$savana$vep_vcf, "_VEP_SAVANA\\.vcf\\.gz$")
  # Severus' VEP pattern must not pick up SAVANA's annotation
  expect_null(inp$callers$severus$vep_vcf)
  expect_null(inp$callers$severus$vcf)

  unlink(file.path(d, "variants/savana/sample.pon_filtered.classified.somatic.vcf"))
  inp = SECTIONS$sv$locate(d, "sample")
  expect_match(inp$callers$savana$vcf, "/sample\\.classified\\.somatic\\.vcf$")
  expect_false(inp$callers$savana$pon_filtered)
})

# ---- Two callers in one table ----------------------------------------------

severus_rec = function(chrom, pos, id, alt, info)
  paste(chrom, pos, id, "N", alt, "60.0", "PASS", info, "GT:VAF:DR:DV", "0/1:0.30:7:3", sep = "\t")

write_gz_vcf = function(path, records) {
  con = gzfile(path, "wt")
  writeLines(c("##fileformat=VCFv4.2",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                       "FORMAT", "SAMPLE1"), collapse = "\t"), records), con)
  close(con)
}

test_that("two callers rbind with a caller column, unique ids and an explicit circos caller", {
  d = stage_sample(c("variants/savana/sample.classified.somatic.vcf" = savana_vcf,
                     "vep/SVs/sample_VEP_SAVANA.vcf.gz"              = savana_vep))
  on.exit(unlink(d, recursive = TRUE))
  sdir = file.path(d, "variants/severus/somatic_SVs")
  dir.create(sdir, recursive = TRUE)
  # A Severus id that collides with a SAVANA one, so the prefixing is exercised
  write_gz_vcf(file.path(sdir, "severus_somatic.vcf.gz"), c(
    severus_rec("chr1", 5000, "severus_DEL1", "<DEL>", "SVTYPE=DEL;END=9000;SVLEN=4000"),
    severus_rec("chr3", 100, "ID_25860_1", "<DUP>", "SVTYPE=DUP;END=900;SVLEN=800")))

  res = SECTIONS$sv$parse(SECTIONS$sv$locate(d, "sample"), list())
  t = res$table
  expect_equal(res$n_by_caller, c(severus = 2L, savana = 12L))
  expect_setequal(unique(t$caller), c("severus", "savana"))
  expect_equal(anyDuplicated(t$id), 0L)
  expect_true(all(startsWith(t[caller == "savana", id], "savana:")))
  expect_true(all(SV_TABLE_COLS %in% names(t)))
  # No Padfoot file anywhere: no Padfoot columns
  expect_false(any(PADFOOT_COLS %in% names(t)))
  expect_false(res$has_padfoot)
  # Severus has records, so it draws the genome-wide circos
  expect_equal(res$circos_caller, "severus")
  expect_equal(nrow(res$circos$nontrans), 2L)
  expect_equal(res$callers$savana$n_records, 20L)
  expect_equal(res$callers$savana$n_typed, 6L)
  # The caller facet appears once two values exist
  expect_match(js_facet_defs(t, cols = "caller"), "savana")
  # Panel matching still runs on the two-caller table
  p = list(has_coords = TRUE, genes = "G", interval_gene = "G",
           chrom = "chr1", start = 635000L, end = 636000L)
  ht = sv_panel_hit_table(t, list(p = p))
  expect_true("savana:ID_25860_1" %in% ht$id)
})

test_that("with SAVANA alone the circos falls back to it", {
  d = stage_sample(c("variants/savana/sample.classified.somatic.vcf" = savana_vcf))
  on.exit(unlink(d, recursive = TRUE))
  res = SECTIONS$sv$parse(SECTIONS$sv$locate(d, "sample"), list())
  expect_equal(res$circos_caller, "savana")
  expect_equal(names(res$n_by_caller), "savana")
  expect_equal(nrow(res$circos$translocations), 3L)
  # Single caller: ids unprefixed
  expect_true("ID_25860_1" %in% res$table$id)
})
