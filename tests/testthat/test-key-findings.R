# Key findings and biomarkers (R/key_findings.R): evidence rules, tiers, panel tagging,
# SV aggregation, the PAR demotion and the biomarker block over absent sections.

kf_fixture = function(f) file.path(repo_root, "tests", "fixtures", "overview", f)

# Two panels as they reach the setup chunk: through the execute_params YAML round trip
kf_panels = function() {
  ps = list(a = load_gene_panel(kf_fixture("panel_a.hg38.tsv"), "hg38"),
            b = load_gene_panel(kf_fixture("panel_b.hg38.tsv"), "hg38"))
  yaml::yaml.load(yaml::as.yaml(ps))
}

kf_vt = function(...) {
  d = data.table(...)
  defaults = list(chrom = "chr17", pos = 7675000L, ref = "C", alt = "T",
                  consequence = "missense_variant", impact = "MODERATE",
                  hgvsp = NA_character_, vaf = 0.4, dp = 30L)
  for (k in names(defaults)) if (!k %in% names(d)) d[, (k) := defaults[[k]]]
  d
}

kf_sv = function(...) {
  d = data.table(...)
  defaults = list(caller = "severus", svtype = "BND", gene_a = NA_character_,
                  gene_b = NA_character_, sv_len = NA_integer_, vaf = 0.3,
                  consequence = NA_character_, impact = NA_character_)
  for (k in names(defaults)) if (!k %in% names(d)) d[, (k) := defaults[[k]]]
  if (!"id" %in% names(d)) d[, id := paste0("sv", seq_len(.N))]
  d
}

kf_run = function(vt = NULL, sd = list(), panels = kf_panels(), ref = "hg38") {
  build_overview(vt, sd, panels, list(purity = NA, ploidy = NA), list(tmb = NA),
                 list(), list(reference = ref, outputs = list(mode = "matched")))$findings
}

# ---- small-variant evidence rules, one at a time ------------------------------------------

test_that("ClinVar pathogenic alone makes a tier-1 finding; Conflicting never does", {
  f = kf_run(kf_vt(symbol = "TP53", impact = "MODIFIER", consequence = "intron_variant",
                   clinvar = "Pathogenic"))
  expect_equal(nrow(f), 1L)
  expect_equal(f$tier, 1L)
  expect_match(f$evidence, "ClinVar: Pathogenic")
  expect_equal(f$panels, "a")

  f = kf_run(kf_vt(symbol = "TP53", impact = "MODIFIER", consequence = "intron_variant",
                   clinvar = "Conflicting_classifications_of_pathogenicity"))
  expect_equal(nrow(f), 0L)
  f = kf_run(kf_vt(symbol = "TP53", impact = "MODIFIER", clinvar = "likely_pathogenic"))
  expect_equal(f$tier, 1L)
})

test_that("HIGH impact alone is tier 1, and shortens the HGVSp", {
  f = kf_run(kf_vt(symbol = "TP53", impact = "HIGH", consequence = "frameshift_variant",
                   hgvsp = "ENSP00000269305.4:p.Pro72fs"))
  expect_equal(f$tier, 1L)
  expect_equal(f$alteration, "p.Pro72fs")
  expect_match(f$evidence, "HIGH impact (frameshift)", fixed = TRUE)
})

test_that("MODERATE needs AlphaMissense likely pathogenic (case-insensitive), at tier 2", {
  f = kf_run(kf_vt(symbol = "TP53", am_class = "Likely_Pathogenic", am_score = 0.98))
  expect_equal(f$tier, 2L)
  expect_match(f$evidence, "AlphaMissense likely pathogenic (0.98)", fixed = TRUE)
  expect_equal(nrow(kf_run(kf_vt(symbol = "TP53", am_class = "likely_benign"))), 0L)
  # No call at all: no finding
  expect_equal(nrow(kf_run(kf_vt(symbol = "TP53"))), 0L)
})

test_that("without an AlphaMissense call, EVE or SIFT+PolyPhen together stand in", {
  f = kf_run(kf_vt(symbol = "TP53", am_class = NA_character_, eve_class = "Pathogenic"))
  expect_equal(nrow(f), 1L)
  expect_match(f$evidence, "EVE pathogenic")
  f = kf_run(kf_vt(symbol = "TP53", sift = "deleterious", polyphen = "probably_damaging"))
  expect_equal(nrow(f), 1L)
  expect_match(f$evidence, "SIFT deleterious + PolyPhen probably damaging", fixed = TRUE)
  # One of the two is not enough
  expect_equal(nrow(kf_run(kf_vt(symbol = "TP53", sift = "deleterious",
                                 polyphen = "possibly_damaging"))), 0L)
  # An AlphaMissense benign call is not overridden by SIFT/PolyPhen
  expect_equal(nrow(kf_run(kf_vt(symbol = "TP53", am_class = "likely_benign", sift = "deleterious",
                                 polyphen = "probably_damaging"))), 0L)
})

test_that("missing optional columns and a NULL table are harmless", {
  expect_equal(nrow(kf_run(NULL)), 0L)
  vt = data.table(symbol = c("TP53", "MYC"), impact = c("HIGH", "LOW"))
  f = kf_run(vt)
  expect_equal(f$gene, "TP53")
  expect_true(is.na(f$vaf))
  expect_equal(nrow(kf_run(data.table(symbol = character(), impact = character()))), 0L)
})

# ---- tiers and roles ------------------------------------------------------------------------

test_that("kf_tier applies role refinements", {
  sv = "small variant"
  expect_equal(kf_tier(sv, "lof", "tsg"), 1L)
  expect_equal(kf_tier(sv, "lof", "oncogene"), 2L)            # LoF in an oncogene
  expect_equal(kf_tier(sv, "clinvar+lof", "oncogene"), 1L)    # ClinVar outranks it
  expect_equal(kf_tier(sv, "missense", "oncogene"), 2L)
  expect_equal(kf_tier(sv, "missense+biallelic", "tsg"), 1L)  # second hit in a TSG
  expect_equal(kf_tier(sv, "missense+biallelic", ""), 2L)     # no role: no adjustment
  expect_equal(kf_tier(sv, "lof", "both"), 1L)
  expect_true(is.na(kf_tier(sv, "", "")))
  expect_equal(kf_tier("CNA", "amplification", "oncogene"), 1L)
  expect_true(is.na(kf_tier("CNA", "amplification", "tsg")))  # amplified TSG: not listed
  expect_equal(kf_tier("CNA", "homdel", "tsg"), 1L)
  expect_equal(kf_tier("CNA", "homdel", "oncogene"), 2L)
  expect_equal(kf_tier("CNA", "homdel+par", "tsg"), 3L)
  expect_equal(kf_tier("SV", "direct", "tsg"), 1L)
  expect_equal(kf_tier("SV", "direct", "oncogene"), 1L)
  expect_equal(kf_tier("SV", "padfoot+padfoot_low", ""), 1L)
  expect_equal(kf_tier("SV", "direct+padfoot_low", ""), 2L)
  expect_equal(kf_tier("SV", "direct+unannotated", ""), 2L)
  expect_equal(kf_tier("SV", "direct+intronic", ""), 3L)
  expect_equal(kf_tier("SV", "direct+recombining", ""), 3L)
  expect_equal(kf_tier("SV", "proximity", ""), 3L)
  # Vectorised, NA role reads as no role
  expect_equal(kf_tier(c(sv, "SV"), c("lof", "proximity"), c(NA, NA)), c(1L, 3L))
  expect_equal(kf_tier_note(sv, "lof", "oncogene"), "loss-of-function in an oncogene")
})

test_that("a truncating variant in a panel oncogene is listed at tier 2 with the reason", {
  f = kf_run(kf_vt(symbol = "MYC", impact = "HIGH", consequence = "stop_gained"))
  expect_equal(f$tier, 2L)
  expect_equal(f$role, "oncogene")
  expect_match(f$evidence, "loss-of-function in an oncogene")
})

test_that("LOH + a listed small variant reads possible biallelic, and lifts a TSG", {
  cna = list(genes = data.table(panel = "a", gene = "TP53", chrom = "chr17", start = 7661779L,
                                end = 7687538L, min_cn = 1, max_cn = 1, min_minor = 0,
                                n_segments = 1L, call = "loss", partial = FALSE, loh = TRUE,
                                source = "wakhan", role = "tsg", expected = 2))
  f = kf_run(kf_vt(symbol = "TP53", am_class = "likely_pathogenic"), sd = list(cna = cna))
  expect_equal(nrow(f), 1L)                 # LOH alone is not a row
  expect_equal(f$tier, 1L)
  expect_match(f$evidence, "possible biallelic (LOH)", fixed = TRUE)
  expect_true(f$loh)
})

# ---- copy number ----------------------------------------------------------------------------

kf_cna = function(gene, chrom, start, end, call, panel = "a", role = NA_character_) {
  list(genes = data.table(panel = panel, gene = gene, chrom = chrom, start = start, end = end,
                          min_cn = if (call == "amplification") 12 else 0,
                          max_cn = if (call == "amplification") 12 else 0, min_minor = 0,
                          n_segments = 1L, call = call, partial = FALSE, loh = FALSE,
                          source = "wakhan", role = role, expected = 2))
}

test_that("a homozygous deletion in a pseudoautosomal region is demoted, on both references", {
  f = kf_run(sd = list(cna = kf_cna("CRLF2", "chrX", 1190490L, 1212723L, "homozygous deletion")))
  expect_equal(f$tier, 3L)
  expect_match(f$evidence, "pseudoautosomal region")
  f = kf_run(sd = list(cna = kf_cna("CRLF2", "chrX", 1000000L, 1020000L, "homozygous deletion")),
             ref = "t2t")
  expect_equal(f$tier, 3L)
  # T2T PAR1 ends at 2,394,410: a gene at 2.5 Mb is outside it there but inside on hg38
  expect_false(kf_in_par("chrX", 2500000, 2510000, "t2t"))
  expect_true(kf_in_par("chrX", 2500000, 2510000, "hg38"))
  expect_true(kf_in_par("chrX", 154000000, 154000100, "t2t"))
  expect_true(kf_in_par("X", 155800000, 155800100, "hg38"))
  expect_false(kf_in_par("chr1", 1000, 2000, "hg38"))
  f = kf_run(sd = list(cna = kf_cna("CDKN2A", "chr9", 21967752L, 21995301L,
                                    "homozygous deletion", panel = "b")))
  expect_equal(f$tier, 1L)
  expect_equal(f$alteration, "homozygous deletion (CN 0)")
})

test_that("an amplified TSG is not listed; an amplified oncogene is tier 1", {
  expect_equal(nrow(kf_run(sd = list(cna = kf_cna("TP53", "chr17", 7661779L, 7687538L,
                                                  "amplification")))), 0L)
  f = kf_run(sd = list(cna = kf_cna("MYC", "chr8", 127735434L, 127742951L, "amplification")))
  expect_equal(f$tier, 1L)
  expect_match(f$alteration, "amplification \\(CN 12\\)")
})

# ---- panel tagging and scopes ---------------------------------------------------------------

test_that("findings are tagged with every panel holding the gene, per_panel_n counts them", {
  sv = kf_sv(svclass = "translocation", chrom_a = "chr22", pos_a = 29280000L,
             chrom_b = "chr12", pos_b = 50800000L)
  o = build_overview(NULL, list(sv = list(table = sv)), kf_panels(), NULL, NULL, list(),
                     list(reference = "hg38"))
  f = o$findings
  # A translocation joining two listed genes is one row; without Padfoot, ordered by chromosome
  expect_equal(nrow(f), 1L)
  expect_equal(f$gene, "ATF1–EWSR1")
  expect_equal(f$genes, "ATF1,EWSR1")
  expect_equal(f$panels, "a,b")
  expect_equal(f$tier, 1L)
  expect_equal(o$per_panel_n, c(a = 1L, b = 1L))
  expect_match(f$alteration, "translocation chr12–chr22")
  expect_false(grepl("partner", f$alteration))
  expect_match(f$evidence, "breakpoints in ATF1 and EWSR1")
})

test_that("the joined row follows Padfoot's gene order, and every other row has genes == gene", {
  sv = kf_sv(svclass = "translocation", chrom_a = "chr22", pos_a = 29280000L,
             chrom_b = "chr12", pos_b = 50800000L, padfoot_impact = "HIGH_ONCO",
             padfoot_genes = "EWSR1 (intron7), ATF1 (intron3)")
  vt = kf_vt(symbol = "TP53", impact = "HIGH", consequence = "stop_gained")
  f = kf_run(vt, sd = list(sv = list(table = sv)))
  expect_setequal(f$gene, c("EWSR1–ATF1", "TP53"))
  expect_equal(f[type == "SV", genes], "EWSR1,ATF1")
  expect_match(f[type == "SV", evidence], "Padfoot HIGH_ONCO: EWSR1 (intron7), ATF1 (intron3)", fixed = TRUE)
  expect_equal(f[type == "small variant", genes], "TP53")
})

test_that("applies_to scopes: an SV-only gene gives no small-variant panel finding and vice versa", {
  # IGH is sv-scoped: a HIGH variant there is not a panel finding (only a non-panel row)
  f = kf_run(kf_vt(symbol = "IGH", impact = "HIGH", consequence = "stop_gained"))
  expect_equal(nrow(f), 1L)
  expect_equal(f$panels, "")
  expect_false(f$in_panel)
  # NOTCH1 is snv-scoped: a breakpoint inside it is not an SV finding
  sv = kf_sv(svclass = "translocation", chrom_a = "chr9", pos_a = 136500000L,
             chrom_b = "chr1", pos_b = 1000000L)
  expect_equal(nrow(kf_run(sd = list(sv = list(table = sv)))), 0L)
  # ...while an IGH breakpoint is
  sv = kf_sv(svclass = "translocation", chrom_a = "chr14", pos_a = 106000000L,
             chrom_b = "chr8", pos_b = 1000000L)
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(f$gene, "IGH")
  expect_equal(f$tier, 1L)
})

test_that("non-panel small variants: only ClinVar/HIGH, capped strongest first", {
  vt = kf_vt(symbol = c("GENE1", "GENE2", "GENE3", "GENE4"),
             impact = c("HIGH", "HIGH", "MODERATE", "HIGH"),
             consequence = c("stop_gained", "stop_gained", "missense_variant", "stop_gained"),
             am_class = c(NA, NA, "likely_pathogenic", NA),
             clinvar = c(NA, "Pathogenic", NA, NA), vaf = c(0.1, 0.2, 0.9, 0.5))
  f = kf_run(vt)
  expect_setequal(f$gene, c("GENE1", "GENE2", "GENE4"))   # predicted missense: panel genes only
  out = kf_small_variant_rows(vt, list(), kf_gene_roles(list()), character(0), max_nonpanel = 2L)
  expect_equal(out$gene, c("GENE2", "GENE4"))              # ClinVar first, then by VAF
  expect_equal(attr(out, "n_nonpanel"), list(kept = 2L, total = 3L))
})

# ---- structural variants --------------------------------------------------------------------

test_that("two callers' calls of one event collapse into one row, listing both callers", {
  sv = kf_sv(id = c("sev_1", "sav_1"), caller = c("severus", "savana"),
             svclass = "translocation", chrom_a = c("chr22", "chr12"),
             pos_a = c(29280000L, 50800003L), chrom_b = c("chr12", "chr22"),
             pos_b = c(50800000L, 29280005L))
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(nrow(f), 1L)
  expect_match(f$alteration, "(SAVANA, Severus)", fixed = TRUE)
  expect_equal(f$source_id, "sev_1,sav_1")
  expect_false(grepl("events", f$evidence))      # one event, however many callers
})

test_that("a Severus intra-chr breakend matching a SAVANA DEL is one deletion row", {
  # CCS15-PB's CDKN2A: the same 2.9 Mb deletion, 1 bp apart, typed only by SAVANA
  sv = kf_sv(id = c("severus_BND1", "ID_1"), caller = c("severus", "savana"),
             svclass = c("intra-chr breakend", "DEL"), svtype = c("BND", "DEL"),
             chrom_a = "chr9", pos_a = c(21969562L, 21969562L),
             chrom_b = "chr9", pos_b = c(24835347L, 24835348L))
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(nrow(f), 1L)
  expect_equal(f$gene, "CDKN2A")
  expect_match(f$alteration, "^deletion chr9:21.97–24.84 Mb \\(SAVANA, Severus\\)$")
  ev = kf_sv_events(sv)
  expect_equal(ev$event, c(1L, 1L))
  expect_equal(ev$ev_class, c("DEL", "DEL"))
})

test_that("calls are separate events when a breakpoint is off, or when one caller made both", {
  base = list(chrom_a = "chr9", chrom_b = "chr9", svclass = "DEL", svtype = "DEL")
  far = do.call(kf_sv, c(base, list(caller = c("severus", "savana"),
                                    pos_a = c(21969562L, 21969562L),
                                    pos_b = c(24835347L, 24835347L + KF_SV_SAME_EVENT_BP + 1L))))
  expect_equal(kf_sv_events(far)$event, c(1L, 2L))
  same_caller = do.call(kf_sv, c(base, list(caller = "severus", pos_a = c(21969562L, 21969570L),
                                            pos_b = c(24835347L, 24835350L))))
  expect_equal(kf_sv_events(same_caller)$event, c(1L, 2L))
  # Two neighbouring 300 bp deletions 400 bp apart are not one event, even across callers
  small = do.call(kf_sv, c(base, list(caller = c("severus", "savana"),
                                      pos_a = c(1000000L, 1000400L), pos_b = c(1000300L, 1000700L))))
  expect_equal(kf_sv_events(small)$event, c(1L, 2L))
})

test_that("a 'near' row from a call matched to a strong event of that gene is dropped", {
  # SAVANA's focal DEL covers CDKN2A; Severus' matching untyped breakend has no breakpoint in
  # it and is only "near" CDKN2A (CCS15-PB's CDKN2B). The two are one event: no tier-3 row
  sv = kf_sv(id = c("sev", "sav"), caller = c("severus", "savana"),
             svclass = c("intra-chr breakend", "DEL"), svtype = c("BND", "DEL"),
             chrom_a = "chr9", pos_a = c(21900000L, 21900001L),
             chrom_b = "chr9", pos_b = c(22100000L, 22100002L))
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(nrow(f), 1L)
  expect_equal(f$tier, 1L)
  expect_match(f$evidence, "focal deletion")
  # ...and the row names both callers of the event, though only SAVANA's call hit it strongly
  expect_match(f$alteration, "(SAVANA, Severus)", fixed = TRUE)
  # Without the SAVANA call the Severus breakend is a tier-3 "near" row on its own
  f = kf_run(sd = list(sv = list(table = sv[1])))
  expect_equal(f$tier, 3L)
})

test_that("a translocation with one listed gene keeps a gene row naming its partner", {
  sv = kf_sv(svclass = "translocation", chrom_a = "chr22", pos_a = 29280000L,
             chrom_b = "chr2", pos_b = 5e6, padfoot_impact = "HIGH_ONCO",
             padfoot_genes = "EWSR1 (intron7), FOO (intron1)")
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(f$gene, "EWSR1")
  expect_equal(f$genes, "EWSR1")
  expect_match(f$alteration, "partner FOO")
})

test_that("rows of one gene are contiguous, groups ordered by their best tier", {
  vt = kf_vt(symbol = c("TP53", "MYC"), impact = c("MODERATE", "HIGH"),
             consequence = c("missense_variant", "stop_gained"),
             am_class = c("likely_pathogenic", NA))
  f = kf_run(vt, sd = list(cna = kf_cna("TP53", "chr17", 7661779L, 7687538L,
                                        "homozygous deletion", role = "tsg")))
  # TP53: tier-1 deletion + tier-2 missense; MYC: tier-2 LoF in an oncogene
  expect_equal(f$gene, c("TP53", "TP53", "MYC"))
  expect_equal(f$tier, c(1L, 2L, 2L))
})

test_that("Padfoot HIGH naming a panel gene is tier 1 even off the panel coordinates", {
  sv = kf_sv(svclass = "translocation", chrom_a = "chr1", pos_a = 5e6, chrom_b = "chr2",
             pos_b = 5e6, padfoot_impact = "HIGH_ONCO",
             padfoot_genes = "EWSR1 (intron7), FOO (intron1)")
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_setequal(f$gene, "EWSR1")
  expect_equal(f$tier, 1L)
  expect_match(f$evidence, "Padfoot HIGH_ONCO")
  expect_match(f$alteration, "partner FOO")
})

test_that("a gene inside a large event is not a finding; a focal deletion covering it is", {
  big = kf_sv(svclass = "INV", svtype = "INV", chrom_a = "chr8", pos_a = 100e6,
              chrom_b = "chr8", pos_b = 140e6)
  expect_equal(nrow(kf_run(sd = list(sv = list(table = big)))), 0L)
  focal = kf_sv(svclass = "DEL", svtype = "DEL", chrom_a = "chr8", pos_a = 127600000L,
                chrom_b = "chr8", pos_b = 127900000L, impact = "HIGH",
                consequence = "transcript_ablation")
  f = kf_run(sd = list(sv = list(table = focal)))
  expect_equal(f$tier, 1L)
  expect_match(f$evidence, "focal deletion")
  # A small intronic deletion annotated by VEP is tier 3; unannotated, tier 2
  intr = kf_sv(svclass = "DEL", svtype = "DEL", chrom_a = "chr17", pos_a = 7670000L,
               chrom_b = "chr17", pos_b = 7670100L, impact = "MODIFIER",
               consequence = "intron_variant")
  expect_equal(kf_run(sd = list(sv = list(table = intr)))$tier, 3L)
  intr[, `:=`(impact = NA_character_, consequence = NA_character_)]
  expect_equal(kf_run(sd = list(sv = list(table = intr)))$tier, 2L)
})

test_that("proximity rows are tier 3, labelled with the distance, and capped per gene", {
  near = 127742951 + 50000   # 50 kb past MYC
  sv = kf_sv(svclass = c("translocation", "DEL", "DUP", "INS", "INV"),
             svtype = c("BND", "DEL", "DUP", "INS", "INV"),
             chrom_a = "chr8", pos_a = near + c(0, 1000, 2000, 3000, 4000),
             chrom_b = c("chr3", "chr8", "chr8", "chr8", "chr8"),
             pos_b = c(1e6, near + 1100, near + 2100, near + 3001, near + 4100))
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_true(all(f$tier == 3L))
  expect_equal(nrow(f), KF_MAX_PROXIMITY_PER_GENE)
  expect_match(f$evidence[1], "near MYC, 50 kb away")
  expect_true(any(grepl("+2 more nearby SV groups not listed", f$evidence, fixed = TRUE)))
})

test_that("a locus-sized deletion inside IGH is V(D)J recombination, tier 3", {
  sv = kf_sv(svclass = "DEL", svtype = "DEL", chrom_a = "chr14", pos_a = 106000000L,
             chrom_b = "chr14", pos_b = 106200000L)
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(f$tier, 3L)
  expect_match(f$evidence, "physiological rearrangement")
})

# ---- biomarkers -----------------------------------------------------------------------------

test_that("every section absent or not run: biomarkers still built, findings empty", {
  sd = list(sv = NULL, cna = NULL, shatterseek = list(status = "not_run"),
            ecdna = list(status = "not_run"), signatures = NULL)
  o = build_overview(NULL, sd, list(), NULL, NULL, NULL, NULL)
  expect_equal(nrow(o$findings), 0L)
  b = o$biomarkers
  expect_setequal(names(b), c("purity", "wgd", "loh", "tmb", "small", "sv", "chromothripsis",
                              "ecdna", "signatures"))
  expect_true(all(vapply(b, function(x) x$status, "") == "not_run"))
  expect_true(all(vapply(b, function(x) x$value, "") == "not run"))
  expect_false(o$tumour_only)
})

test_that("purity names its source, shows the others and flags disagreement", {
  cna = list(source = "wakhan", purity = 0.57, ploidy = 2.1,
             purity_by_source = data.table(source = c("wakhan", "ascat", "savana"),
                                           purity = c(0.57, 1, 0.47), ploidy = c(2.1, 2, 2.2),
                                           used = c(TRUE, FALSE, FALSE)),
             genome = list(wgd = FALSE, loh_fraction = 0.12, wgd_source = "wakhan",
                           major2_fraction = 0.1, ascat_wgd = 1, ascat_loh = 0.2))
  b = build_biomarkers(list(cna = cna), NULL, NULL, list())
  expect_equal(b$purity$value, "0.57 / 2.10")
  expect_match(b$purity$detail, "^Wakhan")
  expect_match(b$purity$detail, "ASCAT 1.00 / 2.00")
  expect_equal(b$purity$status, "warn")
  expect_equal(b$wgd$value, "no")
  expect_match(b$wgd$detail, "ASCAT says yes")
  expect_equal(b$loh$value, "12.0%")
  # No CN section: ASCAT purity/ploidy, named
  b = build_biomarkers(list(), list(purity = 0.8, ploidy = 3.1), NULL, list())
  expect_equal(b$purity$value, "0.80 / 3.10")
  expect_equal(b$purity$detail, "ASCAT")
})

test_that("a whitelist-rescued SV is its own tier, with its read support, below tier 3", {
  sv = kf_sv(id = c("pass", "resc"), svclass = c("DEL", "translocation"),
             svtype = c("DEL", "BND"), chrom_a = c("chr8", "chr22"),
             pos_a = c(127600000L, 29280000L), chrom_b = c("chr8", "chr2"),
             pos_b = c(127900000L, 5e6), impact = c("HIGH", NA),
             consequence = c("transcript_ablation", NA), support = c(12L, 1L),
             filter = c("PASS", "FAIL_MAP_CONS"), wl_rescue = c(NA, "WHITELIST"))
  o = build_overview(NULL, list(sv = list(table = sv)), kf_panels(), NULL, NULL, list(),
                     list(reference = "hg38"))
  f = o$findings
  expect_equal(f$gene, c("MYC", "EWSR1"))
  expect_equal(f$tier, c(1L, KF_TIER_RESCUED))
  expect_equal(f$ev_tier, c(1L, 1L))
  expect_equal(f$rescued, c(FALSE, TRUE))
  expect_match(f$evidence[2], "^whitelist-rescued \\(WHITELIST\\), failed FAIL_MAP_CONS · 1 supporting read · breakpoint in EWSR1")
  # Not counted for the empty-state hints, which describe the default (non-rescued) view
  expect_equal(o$per_panel_n, c(a = 1L, b = 0L))
  # Labelled on the circos only as a rescued row, at its evidence tier
  lab = circos_driver_rows(f)
  expect_equal(lab$label, c("MYC", "EWSR1"))
  expect_equal(lab$tier, c(1L, 1L))
  expect_equal(lab$rescued, c(0L, 1L))
})

test_that("a rescued call that another caller also made is not rescued; tier-3 rescued rows go", {
  # Severus' rescued call and SAVANA's PASS call are one event: it stands with SAVANA's
  sv = kf_sv(id = c("sev", "sav"), caller = c("severus", "savana"),
             svclass = "translocation", chrom_a = c("chr22", "chr22"),
             pos_a = c(29280000L, 29280003L), chrom_b = "chr2", pos_b = c(5e6, 5e6 + 5),
             wl_rescue = c("WHITELIST", NA))
  f = kf_run(sd = list(sv = list(table = sv)))
  expect_equal(nrow(f), 1L)
  expect_false(f$rescued)
  expect_equal(f$tier, 1L)
  # A rescued call only near a gene (tier 3 on its evidence) is not listed at all
  near = kf_sv(svclass = "translocation", chrom_a = "chr8", pos_a = 127742951 + 50000,
               chrom_b = "chr3", pos_b = 1e6, wl_rescue = "WHITELIST")
  expect_equal(nrow(kf_run(sd = list(sv = list(table = near)))), 0L)
  expect_equal(kf_tier("SV", c("direct+rescued", "proximity+rescued", "rescued"), ""),
               c(KF_TIER_RESCUED, KF_TIER_RESCUED, NA_integer_))
})

test_that("the SV biomarker counts whitelist-rescued rows apart", {
  sv = list(n_by_caller = c(severus = 9711L, savana = 57L),
            n_rescued_by_caller = c(severus = 9388L, savana = 0L),
            table = kf_sv(caller = c("severus", "severus", "savana"),
                          svclass = "translocation", chrom_a = "chr1", pos_a = 1L,
                          chrom_b = "chr2", pos_b = 1L, wl_rescue = c("WHITELIST", NA, NA)))
  b = build_biomarkers(list(sv = sv), NULL, NULL, list())
  expect_equal(b$sv$value, "Severus 323 (+9,388 whitelist-rescued) · SAVANA 57")
  expect_match(b$sv$detail, "translocations: Severus 1, SAVANA 1", fixed = TRUE)
})

test_that("SVs are counted per caller, never summed", {
  sv = list(n_by_caller = c(severus = 120L, savana = 64L),
            table = kf_sv(caller = c("severus", "savana", "savana"),
                          svclass = c("translocation", "translocation", "DEL"),
                          chrom_a = "chr1", pos_a = 1L, chrom_b = "chr2", pos_b = 1L))
  b = build_biomarkers(list(sv = sv), NULL, NULL, list())
  expect_equal(b$sv$value, "Severus 120 · SAVANA 64")
  expect_false(grepl("184", b$sv$value))
  # Without n_by_caller, counted from the table
  sv$n_by_caller = NULL
  expect_equal(build_biomarkers(list(sv = sv), NULL, NULL, list())$sv$value,
               "SAVANA 2 · Severus 1")
})

test_that("chromothripsis names chromosome + combination and shows disagreement", {
  ss = list(status = "ok",
            combos = list(severus_ascat = data.table(chrom = "chr1"), savana = data.table(chrom = "chr16")),
            labels = c(severus_ascat = "Severus + ASCAT", savana = "SAVANA"),
            high = data.table(combo = "savana", combo_label = "SAVANA", chrom = "chr16",
                              confidence = "high", criteria = "HC3"))
  b = build_biomarkers(list(shatterseek = ss), NULL, NULL, list())
  expect_equal(b$chromothripsis$value, "chr16")
  expect_match(b$chromothripsis$detail, "chr16 high (SAVANA)", fixed = TRUE)
  expect_match(b$chromothripsis$detail, "disagree")
  expect_equal(b$chromothripsis$status, "warn")
  ss$high = ss$high[0]
  expect_equal(build_biomarkers(list(shatterseek = ss), NULL, NULL, list())$chromothripsis$value,
               "none called")
})

test_that("ecDNA wording follows the classifier, and unclassified is never 'ecDNA'", {
  ec = list(status = "unclassified", amplicons = data.table(amplicon = 1L), n_ecdna = NA, n_bfb = NA)
  b = build_biomarkers(list(ecdna = ec), NULL, NULL, list())
  expect_equal(b$ecdna$value, "focal amplicons found, not classified")
  ec = list(status = "detected", amplicons = data.table(amplicon = 1:2), n_ecdna = 0L, n_bfb = 0L)
  expect_equal(build_biomarkers(list(ecdna = ec), NULL, NULL, list())$ecdna$value, "no ecDNA / BFB")
  ec$n_ecdna = 1L
  expect_equal(build_biomarkers(list(ecdna = ec), NULL, NULL, list())$ecdna$value, "1 ecDNA+")
})

test_that("signatures use only `dominant`, mark artefacts, and say when nothing fitted", {
  sg = list(dominant = data.table(context = c("SBS", "SBS"), signature = c("SBS1", "SBS57"),
                                  mutations = c(100, 50), pct = c(60, 30),
                                  aetiology = c("Clock-like", "Sequencing artefact"),
                                  artefact = c(FALSE, TRUE)))
  b = build_biomarkers(list(signatures = sg), NULL, NULL, list())
  expect_equal(b$signatures$value, "SBS1")
  expect_true(any(grepl("SBS57 30% [artefact]", b$signatures$items, fixed = TRUE)))
  sg$dominant = sg$dominant[0]
  expect_equal(build_biomarkers(list(signatures = sg), NULL, NULL, list())$signatures$value,
               "no well-fitted signature")
})

test_that("the overview reports SBS signatures only, never an indel or doublet fit", {
  sg = list(dominant = data.table(context = c("SBS", "ID", "ID"),
                                  signature = c("SBS5", "ID2", "ID7"),
                                  mutations = c(900, 70, 30), pct = c(100, 70, 30),
                                  aetiology = c("Clock-like", "Replication slippage", "Defective MMR"),
                                  artefact = FALSE))
  b = build_biomarkers(list(signatures = sg), NULL, NULL, list())
  expect_equal(b$signatures$value, "SBS5")
  expect_false(any(grepl("ID7|Defective MMR", b$signatures$items)))
  sg$dominant = sg$dominant[context == "ID"]
  expect_equal(build_biomarkers(list(signatures = sg), NULL, NULL, list())$signatures$value,
               "no well-fitted signature")
})

test_that("tumour-only mode is reported for the caveat", {
  o = build_overview(NULL, list(), list(), NULL, NULL, list(),
                     list(outputs = list(mode = "tumour-only")))
  expect_true(o$tumour_only)
})
