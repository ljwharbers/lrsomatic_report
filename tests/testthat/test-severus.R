# Severus VCF fixtures. Written to a temp .vcf.gz so the real reader (header sniffing,
# INFO parsing, mate collapse) is exercised rather than a hand-built data.table.
write_severus_vcf = function(records, header_extra = character(0)) {
  path = tempfile(fileext = ".vcf.gz")
  con  = gzfile(path, "wt")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##source=Severus_test",
    header_extra,
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "SAMPLE1"), collapse = "\t")
  ), con)
  writeLines(records, con)
  close(con)
  path
}

# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1
rec = function(chrom, pos, id, alt, info, vaf = "0.30") {
  paste(chrom, pos, id, "N", alt, "60.0", "PASS", info, "GT:VAF", paste0("0/1:", vaf),
        sep = "\t")
}

mate_pair = function(prefix = "severus_BND1", chrom_a = "chr1", pos_a = 1000L,
                     chrom_b = "chr2", pos_b = 5000L) {
  c(rec(chrom_a, pos_a, paste0(prefix, "_1"), sprintf("N[%s:%d[", chrom_b, pos_b),
        sprintf("SVTYPE=BND;MATE_ID=%s_2;STRANDS=++", prefix)),
    rec(chrom_b, pos_b, paste0(prefix, "_2"), sprintf("]%s:%d]N", chrom_a, pos_a),
        sprintf("SVTYPE=BND;MATE_ID=%s_1;STRANDS=++", prefix)))
}

test_that("a mate pair collapses to one row carrying both breakend loci", {
  v = write_severus_vcf(mate_pair())
  on.exit(unlink(v))
  r = parse_severus_somatic_records(v)

  expect_equal(nrow(r), 1L)
  expect_equal(r$id, "severus_BND1_1")
  expect_equal(r$id_b, "severus_BND1_2")
  expect_equal(r$chrom_a, "chr1"); expect_equal(r$pos_a, 1000L)
  expect_equal(r$chrom_b, "chr2"); expect_equal(r$pos_b, 5000L)
  expect_equal(r$svclass, "translocation")
})

test_that("an intra-chromosomal breakend is classed apart from a translocation", {
  v = write_severus_vcf(mate_pair(chrom_a = "chr3", pos_a = 100L,
                                  chrom_b = "chr3", pos_b = 900L))
  on.exit(unlink(v))
  expect_equal(parse_severus_somatic_records(v)$svclass, "intra-chr breakend")
})

test_that("an unpaired _1 and an sBND stay as singletons", {
  v = write_severus_vcf(c(
    rec("chr1", 1000L, "severus_BND9_1", "N[chr5:2000[",
        "SVTYPE=BND;MATE_ID=severus_BND9_2"),                  # mate filtered out
    rec("chr7", 7000L, "severus_sBND3", "N.", "SVTYPE=sBND")    # single breakend
  ))
  on.exit(unlink(v))
  r = parse_severus_somatic_records(v)

  expect_equal(nrow(r), 2L)
  expect_true(all(is.na(r$id_b)))
  # The orphaned mate's locus is still named by ALT, so it is not lost.
  expect_equal(r[id == "severus_BND9_1", chrom_b], "chr5")
  expect_equal(r[id == "severus_BND9_1", pos_b], 2000L)
  # A single breakend has no partner in ALT at all.
  expect_true(is.na(r[id == "severus_sBND3", chrom_b]))
  expect_equal(r[id == "severus_sBND3", svclass], "single breakend")
  # sBND must not be filed with the contiguous types, whose track has no place for it.
  expect_equal(nrow(severus_circos_tracks(r)$nontrans), 0L)
})

test_that("the partner locus parses from all four ALT bracket forms, with or without chr", {
  alts = c("N[chr7:24547089[", "]chr7:24547089]N", "N]chr7:24547089]", "[chr7:24547089[N")
  for (a in alts) {
    v = write_severus_vcf(rec("chr1", 100L, "b1", a, "SVTYPE=BND"))
    r = parse_severus_somatic_records(v)
    expect_equal(r$chrom_b, "chr7", info = a)
    expect_equal(r$pos_b, 24547089L, info = a)
    unlink(v)
  }
  # A bare-contig reference: the "chr" prefix is added, not required.
  v = write_severus_vcf(rec("1", 100L, "b1", "N[7:555[", "SVTYPE=BND"))
  on.exit(unlink(v))
  r = parse_severus_somatic_records(v)
  expect_equal(r$chrom_a, "chr1")
  expect_equal(r$chrom_b, "chr7")
  expect_equal(r$pos_b, 555L)
})

test_that("a VCF with no END and a malformed BND ALT parses instead of erroring", {
  # No record carries END= or SVLEN=, and one ALT has no partner locus at all.
  v = write_severus_vcf(c(
    rec("chr1", 100L, "i1", "ACGT", "SVTYPE=INS"),
    rec("chr2", 200L, "b1", "N.",   "SVTYPE=BND")
  ))
  on.exit(unlink(v))
  r = parse_severus_somatic_records(v)

  expect_equal(nrow(r), 2L)
  expect_equal(r[id == "i1", pos_b], 100L)   # END falls back to POS
  expect_true(is.na(r[id == "b1", pos_b]))
  expect_true(all(is.na(r$sv_len)))
})

test_that("VAF comes from the FORMAT field", {
  v = write_severus_vcf(rec("chr1", 100L, "d1", "<DEL>", "SVTYPE=DEL;END=900;SVLEN=800",
                            vaf = "0.42"))
  on.exit(unlink(v))
  r = parse_severus_somatic_records(v)
  expect_equal(r$vaf, 0.42)
  expect_equal(r$pos_b, 900L)
  expect_equal(r$sv_len, 800L)
})

test_that("a VCF with no FORMAT/SAMPLE columns still yields its SVs", {
  # Reading a fixed 10 columns would fail here and return zero rows, which reads
  # exactly like "this sample has no SVs".
  path = tempfile(fileext = ".vcf.gz")
  con  = gzfile(path, "wt")
  writeLines(c("##fileformat=VCFv4.2",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
                     collapse = "\t"),
               paste("chr1", 100, "d1", "N", "<DEL>", "60.0", "PASS",
                     "SVTYPE=DEL;END=900", sep = "\t")), con)
  close(con)
  on.exit(unlink(path))

  r = parse_severus_somatic_records(path)
  expect_equal(nrow(r), 1L)
  expect_equal(r$pos_b, 900L)
  expect_true(is.na(r$vaf))
})

test_that("circos links are one per rearrangement", {
  v = write_severus_vcf(c(mate_pair("severus_BNDa"),
                          mate_pair("severus_BNDb", chrom_b = "chr4", pos_b = 8000L)))
  on.exit(unlink(v))
  tr = severus_circos_tracks(parse_severus_somatic_records(v))
  expect_equal(nrow(tr$translocations), 2L)
})

# ---- VEP annotation join -------------------------------------------------

write_vep_sv_vcf = function(records) {
  path = tempfile(fileext = ".vcf.gz")
  con  = gzfile(path, "wt")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele|Consequence|IMPACT|SYMBOL|Gene\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), collapse = "\t")
  ), con)
  writeLines(records, con)
  close(con)
  path
}

vep_rec = function(chrom, pos, id, alt, csq) {
  paste(chrom, pos, id, "N", alt, "60.0", "PASS", paste0("SVTYPE=X;CSQ=", csq), sep = "\t")
}

test_that("the annotation join is keyed on record ID, not locus", {
  # A DEL and a BND start at the same base: keying on the locus leaks the DEL's gene
  # onto the breakend (and vice versa).
  sv = write_severus_vcf(c(
    rec("chr1", 5000L, "severus_DEL1", "<DEL>", "SVTYPE=DEL;END=9000;SVLEN=4000"),
    rec("chr1", 5000L, "severus_BND2_1", "N[chr9:700[", "SVTYPE=BND;MATE_ID=severus_BND2_2"),
    rec("chr9", 700L,  "severus_BND2_2", "]chr1:5000]N", "SVTYPE=BND;MATE_ID=severus_BND2_1")
  ))
  vep = write_vep_sv_vcf(c(
    vep_rec("chr1", 5000L, "severus_DEL1",   "<DEL>",       "deletion|feature_truncation|HIGH|DELGENE|ENSG1"),
    vep_rec("chr1", 5000L, "severus_BND2_1", "N[chr9:700[", "N|feature_truncation|MODIFIER|SIDEA|ENSG2"),
    vep_rec("chr9", 700L,  "severus_BND2_2", "]chr1:5000]N", "N|feature_truncation|HIGH|SIDEB|ENSG3")
  ))
  on.exit(unlink(c(sv, vep)))

  t = build_sv_table_from_vep(sv, vep)
  expect_equal(nrow(t), 2L)
  expect_equal(t[id == "severus_DEL1", gene_a], "DELGENE")
  expect_equal(t[id == "severus_BND2_1", gene_a], "SIDEA")
  # Side B's symbols come from the mate record, which is what makes a panel gene at the
  # far breakend visible at all.
  expect_equal(t[id == "severus_BND2_1", gene_b], "SIDEB")
  # ...and the higher-impact side supplies consequence/impact.
  expect_equal(t[id == "severus_BND2_1", impact], "HIGH")
})

test_that("a missing SV VEP VCF leaves the gene columns empty but keeps every SV", {
  sv = write_severus_vcf(mate_pair())
  on.exit(unlink(sv))
  t = build_sv_table_from_vep(sv, NULL)
  expect_equal(nrow(t), 1L)
  expect_true(is.na(t$gene_a) && is.na(t$gene_b))
  expect_equal(t$pos_b, 5000L)
})

# ---- Panel matching ------------------------------------------------------

sv_rows = function() data.table(
  id      = c("bnd1", "del1"),
  svtype  = c("BND", "DEL"),
  svclass = c("translocation", "DEL"),
  chrom_a = c("chr8", "chr10"),
  pos_a   = c(1000000L, 5000000L),
  chrom_b = c("chr14", "chr10"),
  pos_b   = c(2000000L, 5100000L),
  gene_a  = c("SOMEGENE", NA_character_),
  gene_b  = c(NA_character_, NA_character_)
)

coord_panel = function(gene, chrom, start, end)
  list(has_coords = TRUE, genes = gene, interval_gene = gene,
       chrom = chrom, start = as.integer(start), end = as.integer(end))

test_that("a coordinate panel matches exactly at the window edge and not beyond it", {
  sv = sv_rows()
  # BND side A at chr8:1,000,000 — an interval starting exactly 1 Mb away matches.
  expect_equal(sv_panel_hits(sv, coord_panel("EDGE", "chr8", 2000000, 2000100))[1],
               "EDGE (A)")
  expect_equal(sv_panel_hits(sv, coord_panel("PAST", "chr8", 2000001, 2000100))[1], "")
  # DEL span chr10:5,000,000-5,100,000 — 100 kb window on the span, not on each end.
  expect_equal(sv_panel_hits(sv, coord_panel("EDGE", "chr10", 5200000, 5200100))[2],
               "EDGE (span)")
  expect_equal(sv_panel_hits(sv, coord_panel("PAST", "chr10", 5200001, 5200100))[2], "")
})

test_that("a panel interval matching only side B is still a hit", {
  sv = sv_rows()
  hits = sv_panel_hits(sv, coord_panel("FARSIDE", "chr14", 2000500, 2001000))
  expect_equal(hits[1], "FARSIDE (B)")
  expect_equal(hits[2], "")
})

test_that("a symbol-only panel matches either side's annotated gene, with no window", {
  sv = sv_rows()
  sym = list(has_coords = FALSE, genes = c("SOMEGENE", "OTHER"))
  expect_equal(sv_panel_hits(sv, sym), c("SOMEGENE (A)", ""))
  # A gene 1 Mb away is invisible to symbol matching — the failure this replaces.
  expect_equal(sv_panel_hits(sv, list(has_coords = FALSE, genes = "EDGE")), c("", ""))
})

test_that("no panel means no hits, and an empty table is handled", {
  expect_equal(sv_panel_hits(sv_rows(), NULL), c("", ""))
  expect_equal(sv_panel_hits(data.table(), coord_panel("X", "chr1", 1, 2)), character(0))
})

test_that("build_sv_table maps the gene-annotated TSV onto the shared schema", {
  tsv = data.table(
    ID = c("SV1", "SV2"), SVTYPE = c("DEL", "BND"), DETAILED_TYPE = c("DEL", "BND"),
    START_CHROM = c("1", "chr2"), START_POS = c(1000L, 2000L),
    END_CHROM   = c("1", "chr7"), END_POS   = c(5000L, 6000L),
    SV_LEN = c(4000L, NA_integer_), VAF = c(0.4, 0.5),
    NHL_GENE_HITS = c("MYC;BCL2", NA_character_)
  )
  t = build_sv_table(tsv)
  expect_equal(nrow(t), 2L)
  expect_equal(t$chrom_a, c("chr1", "chr2"))   # prefix normalised
  expect_equal(t$svclass, c("DEL", "translocation"))
  expect_equal(t$gene_a[1], "MYC;BCL2")
  # Coordinate matching works on this path too, because it has both loci.
  expect_equal(sv_panel_hits(t, coord_panel("NEAR", "chr7", 6000, 6100))[2], "NEAR (B)")
})
