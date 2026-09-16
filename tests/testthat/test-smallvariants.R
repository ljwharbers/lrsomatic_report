# Minimal synthetic VEP data (mimics parse_vep_text output)
make_vep = function() {
  data.table(
    chrom        = c("chr1", "chr1", "chr2"),
    pos          = c(100L, 200L, 300L),
    ref          = c("A", "C", "G"),
    alt          = c("T", "G", "A"),
    symbol       = c("MYC", "TP53", "KRAS"),
    gene_id      = c("ENSG001", "ENSG002", "ENSG003"),
    consequence  = c("missense_variant", "stop_gained", "synonymous_variant"),
    impact       = c("MODERATE", "HIGH", "LOW"),
    hgvsp        = c("p.A1T", "p.Q2*", NA_character_),
    existing     = c(NA_character_, NA_character_, NA_character_),
    dbsnp        = c(NA_character_, NA_character_, NA_character_),
    cosmic       = c(NA_character_, NA_character_, NA_character_),
    sift         = c(NA_character_, NA_character_, NA_character_),
    polyphen     = c(NA_character_, NA_character_, NA_character_),
    # parse_vep_text() reads coordinates from `variant_id`, i.e. VEP's own space.
    coord_space  = "vep"
  )
}

test_that("build_variant_table with NULL panel returns all variants", {
  vep = make_vep()
  result = build_variant_table(vep, NULL, gene_panel = NULL)
  expect_equal(nrow(result), 3L)
  expect_true(all(c("MYC", "TP53", "KRAS") %in% result$symbol))
})

test_that("build_variant_table with panel filters to panel genes", {
  vep = make_vep()
  result = build_variant_table(vep, NULL, gene_panel = c("MYC", "TP53"))
  expect_equal(nrow(result), 2L)
  expect_false("KRAS" %in% result$symbol)
})

test_that("build_variant_table with empty panel returns empty table", {
  vep = make_vep()
  result = build_variant_table(vep, NULL, gene_panel = character(0))
  expect_equal(nrow(result), 0L)
})

# --- parse_caller_vcf(): FORMAT tag and sample column -----------------------

# A gzipped VCF built from explicit FORMAT/sample strings, so each test states exactly
# the layout it is about.
write_vcf = function(fmt, samples, records, extra_header = character(0)) {
  path = tempfile(fileext = ".vcf.gz")
  head = paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                 "FORMAT", samples), collapse = "\t")
  body = vapply(records, function(r) paste(c(r$fixed, fmt, r$samples), collapse = "\t"),
                character(1))
  con = gzfile(path, "wb")
  writeLines(c("##fileformat=VCFv4.2", extra_header, head, body), con)
  close(con)
  path
}

snv_rec = function(...) list(fixed = c("chr1", "100", ".", "A", "G", ".", "PASS", "."),
                             samples = c(...))

test_that("parse_caller_vcf reads the allele fraction from AF or VAF", {
  # LRSomatic renames this tag by prioritised caller (STANDARDIZE_AF), so both reach us.
  af  = parse_caller_vcf(write_vcf("GT:DP:AF",  "TUM", list(snv_rec("0/1:50:0.42"))))
  vaf = parse_caller_vcf(write_vcf("GT:DP:VAF", "TUM", list(snv_rec("0/1:50:0.42"))))
  expect_equal(af$vaf, 0.42)
  expect_equal(vaf$vaf, 0.42)
})

test_that("parse_caller_vcf derives the allele fraction from AD when neither tag exists", {
  res = parse_caller_vcf(write_vcf("GT:DP:AD", "TUM", list(snv_rec("0/1:40:30,10"))))
  expect_equal(res$vaf, 0.25)
})

test_that("parse_caller_vcf reads the tumour sample, not whichever sorts first", {
  # ClairS matched mode emits two samples; reading column 10 positionally would report
  # the normal's VAF for every somatic variant, with no error to reveal it.
  f = write_vcf("GT:DP:VAF", c("NORMAL", "TUMOR"),
                list(snv_rec("0/0:60:0.01", "0/1:50:0.44")))
  expect_equal(parse_caller_vcf(f, sample_id = "TUMOR")$vaf, 0.44)
  expect_equal(parse_caller_vcf(f, sample_id = "TUMOR")$dp, 50L)
  # With no sample_id to match, fall back to the sample that is not named like a normal.
  expect_equal(suppressWarnings(parse_caller_vcf(f))$vaf, 0.44)
})

test_that("parse_caller_vcf returns NULL for a sites-only VCF instead of erroring", {
  path = tempfile(fileext = ".vcf.gz")
  con = gzfile(path, "wb")
  writeLines(c("##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               "chr1\t100\t.\tA\tG\t.\tPASS\t."), con)
  close(con)
  expect_null(parse_caller_vcf(path))
})

# --- coordinate space: the CSQ path must join its indels ---------------------

# A CSQ-format (VEP --vcf) file and the caller VCF it was annotated from, carrying the
# same three variants in the same VCF representation.
csq_variants = list(
  list(fixed = c("chr1", "100", ".", "A",  "G",  ".", "PASS"), csq = "G|missense_variant|MODERATE|TP53|E1"),
  list(fixed = c("chr1", "200", ".", "TG", "T",  ".", "PASS"), csq = "-|frameshift_variant|HIGH|BRCA1|E2"),
  list(fixed = c("chr1", "300", ".", "A",  "AT", ".", "PASS"), csq = "T|frameshift_variant|HIGH|MYC|E3"))

write_gz = function(lines) {
  path = tempfile(fileext = ".vcf.gz")
  con = gzfile(path, "wb"); writeLines(lines, con); close(con)
  path
}

test_that("indels on the CSQ path join their VAF (coordinate space is declared, not assumed)", {
  head = paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                 "FORMAT", "TUM"), collapse = "\t")
  vep_f = write_gz(c(
    "##fileformat=VCFv4.2",
    paste0("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: ",
           "Allele|Consequence|IMPACT|SYMBOL|Gene\">"),
    head,
    vapply(csq_variants, function(v)
      paste(c(v$fixed, paste0("CSQ=", v$csq), "GT", "0/1"), collapse = "\t"), character(1))))
  vaf_f = write_gz(c(
    "##fileformat=VCFv4.2", head,
    vapply(csq_variants, function(v)
      paste(c(v$fixed, ".", "GT:DP:VAF", "0/1:50:0.33"), collapse = "\t"), character(1))))

  vep = parse_vep_vcf(vep_f)
  expect_equal(unique(vep$coord_space), "vcf")   # VCF POS, not VEP's shifted one

  vaf = parse_caller_vcf(vaf_f, "somatic", sample_id = "TUM")
  res = build_variant_table(vep, vaf)

  # All three classes, not just the SNV: keying the VCF-space VEP table as VEP-space
  # shifted one side of the key and no indel could ever match.
  expect_equal(sum(!is.na(res$vaf)), 3L)
  expect_true(all(res[nchar(ref) != nchar(alt)]$vaf == 0.33))
})

test_that("the VEP text path still joins, and declares VEP space", {
  vep = make_vep()
  expect_equal(unique(vep$coord_space), "vep")
  # VEP text reports an indel one base right of the VCF anchor, in dash form.
  vep = rbind(vep, data.table(
    chrom = "chr3", pos = 501L, ref = "G", alt = "-", symbol = "RB1",
    gene_id = "ENSG004", consequence = "frameshift_variant", impact = "HIGH",
    hgvsp = NA_character_, existing = NA_character_, dbsnp = NA_character_,
    cosmic = NA_character_, sift = NA_character_, polyphen = NA_character_,
    coord_space = "vep"), fill = TRUE)

  vaf = data.table(chrom = c("chr1", "chr3"), pos = c(100L, 500L),
                   ref = c("A", "CG"), alt = c("T", "C"),
                   vaf = c(0.4, 0.2), dp = 50L, gt = "0/1", ps = NA_character_)
  res = build_variant_table(vep, vaf)
  expect_equal(res[symbol == "MYC"]$vaf, 0.4)
  expect_equal(res[symbol == "RB1"]$vaf, 0.2)   # deletion, reconciled by variant_key()
})

# --- VEP plugin annotations --------------------------------------------------

# A CSQ VCF with the given Format fields, one entry per record, plus optional header lines
write_csq_vcf = function(format, entries, extra_header = character(0)) {
  head = paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), collapse = "\t")
  write_gz(c(
    "##fileformat=VCFv4.2", extra_header,
    paste0("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: ", format, "\">"),
    head,
    vapply(seq_along(entries), function(i)
      paste(c("chr1", 100L * i, ".", "A", "G", ".", "PASS", paste0("CSQ=", entries[i])),
            collapse = "\t"), character(1))))
}
csq_base = "Allele|Consequence|IMPACT|SYMBOL|Gene"

test_that("the PolyPhen_SIFT, AlphaMissenseProtein and ClinVar plugin fields (T2T) become columns", {
  f = write_csq_vcf(
    paste(csq_base, "SIFT_pred|SIFT_score|PolyPhen_humvar_pred|PolyPhen_humvar_score",
          "AlphaMissenseProtein_match|am_class|am_pathogenicity|ClinVar|ClinVar_CLNSIG", sep = "|"),
    c("G|missense_variant|MODERATE|TP53|E1|deleterious|0|probably_damaging|0.999|gene_aa|pathogenic|0.97|12345|Pathogenic&risk_factor",
      "G|missense_variant|MODERATE|KRAS|E2|||||aa_mismatch||||"),
    extra_header = '##VEP="v115.2" API="v115" assembly="T2T-CHM13v2.0"')
  v = parse_vep_vcf(f)
  expect_equal(v$sift,           c("deleterious", NA));       expect_equal(v$sift_score,     c(0, NA))
  expect_equal(v$polyphen,       c("probably_damaging", NA)); expect_equal(v$polyphen_score, c(0.999, NA))
  expect_equal(v$am_class,       c("pathogenic", NA));        expect_equal(v$am_score,       c(0.97, NA))
  expect_equal(v$clinvar,        c("Pathogenic,risk_factor", NA))   # "&" rewritten so the facet can split it
  expect_equal(v$clinvar_id,     c("12345", NA))

  # Plugins this run did not have get no column at all, and the T2T cache carries no dbSNP/COSMIC
  expect_false(any(c("cadd", "revel", "eve_class", "eve_score", "dbsnp", "cosmic") %in% names(v)))
  vf = attr(v, "vep_fields")
  expect_equal(unname(vf["sift"]), "SIFT_pred")
  expect_true(is.na(vf["cadd"])); expect_true(is.na(vf["dbsnp"]))

  # The display table carries the columns and the attribute through merge()/unique()
  tb = build_variant_table(v, NULL)
  expect_true(all(c("sift", "sift_score", "am_class", "am_score", "clinvar", "clinvar_id") %in% names(tb)))
  expect_false("cadd" %in% names(tb))
  expect_equal(attr(tb, "vep_fields"), vf)
})

test_that("built-in SIFT/PolyPhen `pred(score)` splits into the same class and score columns (hg38)", {
  f = write_csq_vcf(
    paste(csq_base, "SIFT|PolyPhen|CLIN_SIG|am_class|am_pathogenicity|CADD_PHRED|CADD_RAW|REVEL|EVE_CLASS|EVE_SCORE", sep = "|"),
    c("G|missense_variant|MODERATE|TP53|E1|deleterious(0.01)|benign(0)|pathogenic|likely_pathogenic|0.98|25.3|4.1|0.9|Pathogenic|0.7375145263264659",
      "G|missense_variant|MODERATE|KRAS|E2|deleterious_low_confidence(0)|||||1.2|-0.1|||"),
    extra_header = '##VEP="v115.2" API="v115" assembly="GRCh38.p14" COSMIC="99" dbSNP="156"')
  v = parse_vep_vcf(f)
  expect_equal(v$sift,       c("deleterious", "deleterious_low_confidence"))
  expect_equal(v$sift_score, c(0.01, 0))
  expect_equal(v$polyphen,   c("benign", NA)); expect_equal(v$polyphen_score, c(0, NA))
  expect_equal(v$clinvar,    c("pathogenic", NA))   # the cache's CLIN_SIG when no ClinVar --custom track ran
  expect_false("clinvar_id" %in% names(v))
  expect_equal(v$cadd,       c(25.3, 1.2)); expect_equal(v$revel, c(0.9, NA))
  expect_equal(v$eve_class,  c("Pathogenic", NA)); expect_equal(v$eve_score[1], 0.7375145263264659)
  expect_true(all(c("dbsnp", "cosmic") %in% names(v)))
  vf = attr(v, "vep_fields")
  expect_equal(unname(vf[c("sift", "sift_score", "clinvar", "dbsnp")]),
               c("SIFT", "SIFT", "CLIN_SIG", "Existing_variation"))
})

test_that("a CSQ header with no prediction field yields no prediction column; an undeclared cache keeps dbSNP/COSMIC", {
  v = parse_vep_vcf(write_csq_vcf(csq_base, "G|missense_variant|MODERATE|TP53|E1"))
  expect_false(any(names(VEP_ANNOTATION_FIELDS) %in% names(v)))
  expect_true(all(c("dbsnp", "cosmic") %in% names(v)))   # no ##VEP line, so nothing says the cache lacks them
  s = vep_annotation_summary(attr(v, "vep_fields"))
  expect_equal(s[present == TRUE]$source, c("dbSNP", "COSMIC"))
  expect_equal(s[source == "SIFT"]$fields, "")
})

test_that("the text path resolves the same plugin keys from Extra, whether declared in the header or only present in a row", {
  cols = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
           "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
           "Codons", "Existing_variation", "Extra")
  row = function(vid, loc, allele, gene, csq, extra)
    paste(c(vid, loc, allele, gene, "T1", "Transcript", csq, "-", "-", "-", "-", "-", "-", extra),
          collapse = "\t")
  f = write_gz(c(
    "## ENSEMBL VARIANT EFFECT PREDICTOR v115.2",
    "## Output produced at 2026-09-15 20:41:02",
    "## Extra column keys:",
    "## IMPACT : Subtype of consequence type",
    "## SYMBOL : Gene symbol",
    "## SIFT_pred : SIFT prediction",
    "## SIFT_score : SIFT score",
    paste0("#", paste(cols, collapse = "\t")),
    row("chr1_100_A/G", "chr1:100", "G", "E1", "missense_variant",
        "IMPACT=MODERATE;SYMBOL=TP53;SIFT_pred=deleterious;SIFT_score=0.02;am_class=pathogenic;am_pathogenicity=0.9"),
    row("chr1_200_C/T", "chr1:200", "T", "E2", "synonymous_variant", "IMPACT=LOW;SYMBOL=KRAS")))
  v = parse_vep_text(f)
  expect_equal(v$symbol,   c("TP53", "KRAS"))
  expect_equal(v$sift,     c("deleterious", NA)); expect_equal(v$sift_score, c(0.02, NA))
  expect_equal(v$am_class, c("pathogenic", NA));  expect_equal(v$am_score,   c(0.9, NA))
  expect_false("polyphen" %in% names(v))
  # A text header names its dbSNP/COSMIC releases; this one has neither line
  expect_false(any(c("dbsnp", "cosmic") %in% names(v)))
})

test_that("the text path reads Existing_variation from its own column, not from Extra", {
  cols = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
           "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
           "Codons", "Existing_variation", "Extra")
  row = function(vid, loc, allele, existing)
    paste(c(vid, loc, allele, "E1", "T1", "Transcript", "missense_variant",
            "-", "-", "-", "-", "-", existing, "IMPACT=MODERATE;SYMBOL=TP53"), collapse = "\t")
  f = write_gz(c(
    "## ENSEMBL VARIANT EFFECT PREDICTOR v115.2",
    "## COSMIC version 99",
    "## dbSNP version 156",
    "## Extra column keys:",
    "## IMPACT : Subtype of consequence type",
    "## SYMBOL : Gene symbol",
    "## CLIN_SIG : ClinVar clinical significance of the dbSNP variant",
    paste0("#", paste(cols, collapse = "\t")),
    row("chr1_100_A/G", "chr1:100", "G", "rs886974256"),
    row("chr1_200_C/T", "chr1:200", "T", "rs772325487,COSV100527437"),
    row("chr1_300_G/A", "chr1:300", "A", "-")))
  v = parse_vep_text(f)
  expect_equal(v$existing, c("rs886974256", "rs772325487,COSV100527437", NA))
  expect_equal(v$dbsnp,    c("rs886974256", "rs772325487", NA))
  expect_equal(v$cosmic,   c(NA, "COSV100527437", NA))
})

test_that("vep_cache_sources reads both header shapes and stays agnostic without a VEP line", {
  expect_equal(vep_cache_sources('##VEP="v115" dbSNP="156" COSMIC="99"'),   c(dbsnp = TRUE,  cosmic = TRUE))
  expect_equal(vep_cache_sources('##VEP="v115" assembly="T2T-CHM13v2.0"'), c(dbsnp = FALSE, cosmic = FALSE))
  expect_equal(vep_cache_sources(c("## ENSEMBL VARIANT EFFECT PREDICTOR v115", "## dbSNP 156")),
               c(dbsnp = TRUE, cosmic = FALSE))
  expect_equal(vep_cache_sources(c("## ENSEMBL VARIANT EFFECT PREDICTOR v115",
                                   "## dbSNP version 156", "## COSMIC version 99")),
               c(dbsnp = TRUE, cosmic = TRUE))
  expect_equal(vep_cache_sources("##fileformat=VCFv4.2"), c(dbsnp = NA, cosmic = NA))
})

test_that("a VEP text header describing its Extra keys does not read as a dbSNP-carrying cache", {
  # Every VEP text header carries this line, whether or not the cache has dbSNP at all
  t2t = c("## ENSEMBL VARIANT EFFECT PREDICTOR v115.2",
          "## assembly version T2T-CHM13v2.0",
          "## Extra column keys:",
          "## CLIN_SIG : ClinVar clinical significance of the dbSNP variant",
          "## SOMATIC : Somatic status of existing variant")
  expect_equal(vep_cache_sources(t2t), c(dbsnp = FALSE, cosmic = FALSE))
})

test_that("extract_extra_keys matches extract_extra_key and keeps '=' inside a value", {
  ex = c("SYMBOL=TP53;HGVSp=ENSP1:p.A1T;NOTE=a=b", "-", NA, "SYMBOL=KRAS")
  out = extract_extra_keys(ex, c("SYMBOL", "NOTE", "MISSING"))
  expect_equal(out$SYMBOL,  c("TP53", NA, NA, "KRAS"))
  expect_equal(out$NOTE,    c("a=b", NA, NA, NA))
  expect_equal(out$MISSING, rep(NA_character_, 4))
  expect_equal(extract_extra_key(ex, "SYMBOL"), out$SYMBOL)
})

# --- vaf_provenance() -------------------------------------------------------

# Write a minimal gzipped VCF header (plus one record, so the scan has to stop on its own)
make_vcf_gz = function(path, source_line = NULL) {
  lines = c("##fileformat=VCFv4.2")
  if (!is.null(source_line)) lines = c(lines, paste0("##source=", source_line))
  lines = c(lines,
            "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
            "chr1\t100\t.\tA\tT\t30\tPASS\t.\tGT:AF\t0/1:0.4")
  con = gzfile(path, "wb")
  writeLines(lines, con)
  close(con)
  path
}

test_that("vaf_provenance returns NULL when there is nothing to describe", {
  expect_null(vaf_provenance(NULL))
  expect_null(vaf_provenance(character(0)))
  expect_null(vaf_provenance(NA_character_))
  expect_null(vaf_provenance(""))
})

test_that("vaf_provenance reports the path and the declared ##source", {
  d = withr::local_tempdir()
  f = make_vcf_gz(file.path(d, "somatic.vcf.gz"), source_line = "Clair3")

  prov = vaf_provenance(f, sample_dir = d)
  expect_equal(prov$paths, "somatic.vcf.gz")   # relative to sample_dir
  expect_equal(prov$sources, "Clair3")
})

test_that("vaf_provenance keeps the absolute path when it is outside sample_dir", {
  d = withr::local_tempdir()
  f = make_vcf_gz(file.path(d, "somatic.vcf.gz"), source_line = "Clair3")

  prov = vaf_provenance(f, sample_dir = file.path(d, "elsewhere"))
  expect_equal(prov$paths, f)
})

test_that("vaf_provenance handles several files and a missing ##source", {
  d = withr::local_tempdir()
  snvs  = make_vcf_gz(file.path(d, "snvs.vcf.gz"),  source_line = "ClairS")
  indel = make_vcf_gz(file.path(d, "indel.vcf.gz"), source_line = NULL)

  prov = vaf_provenance(c(snvs, indel), sample_dir = d)
  expect_equal(prov$paths, c("snvs.vcf.gz", "indel.vcf.gz"))
  expect_equal(prov$sources, "ClairS")   # deduped; the second file declares none
})

test_that("vaf_provenance does not fail on a path that does not exist", {
  prov = vaf_provenance("/nonexistent/somatic.vcf.gz")
  expect_equal(prov$paths, "/nonexistent/somatic.vcf.gz")
  expect_length(prov$sources, 0L)
})
