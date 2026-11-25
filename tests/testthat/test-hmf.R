# HMF driver / fusion / hotspot annotation (R/hmf.R). Fixtures in tests/fixtures/hmf/ are
# synthetic (real gene names, made-up flags and coordinates): see their README
hmf_fx = file.path(repo_root, "tests", "fixtures", "hmf")
hmf_files = function(...) {
  f = list(driver_panel = file.path(hmf_fx, "DriverGenePanel.38.tsv"),
           ensembl_data_dir = file.path(hmf_fx, "ensembl_data"),
           known_fusions = file.path(hmf_fx, "known_fusion_data.38.csv"),
           hotspots = file.path(hmf_fx, "KnownHotspots.somatic.38.vcf"),
           actionable_bed = file.path(hmf_fx, "ActionableCodingPanel.38.bed"))
  modifyList(f, list(...))
}
hmf_res = function(...) hmf_load_resources(hmf_files(...))

# A v3-style copy of a fixture: v3 renames IG_* to ENHANCER_* and adds panel columns
hmf_v3_copy = function(name, edit) {
  out = tempfile(fileext = paste0(".", tools::file_ext(name)))
  writeLines(edit(readLines(file.path(hmf_fx, name))), out)
  out
}

# ---- Loaders ---------------------------------------------------------------------------

test_that("the driver panel reads booleans, roles and a defaulted amplification ratio", {
  p = hmf_read_panel(file.path(hmf_fx, "DriverGenePanel.38.tsv"))
  expect_equal(nrow(p), 7L)
  expect_equal(attr(p, "n_columns"), 15L)
  expect_true(p[gene == "BRAF", report_missense & report_hotspot & report_amplification])
  expect_false(p[gene == "BRAF", report_nonsense])
  expect_equal(p[gene == "TP53", role], "tsg")
  expect_equal(p[gene == "MYC", role], "oncogene")
  expect_true(all(p$amp_ratio == HMF_AMP_RATIO_DEFAULT))
})

test_that("a v3 panel's extra columns are read, and absent flags read as false", {
  v3 = hmf_v3_copy("DriverGenePanel.38.tsv", function(l) {
    l[1] = paste0(l[1], "\tamplificationRatio\treportLoh")
    l[-1] = paste0(l[-1], ifelse(startsWith(l[-1], "MYC\t"), "\t2\tfalse", "\t\tfalse"))
    # drop reportSplice: a targeted panel may lack a column
    vapply(strsplit(l, "\t", fixed = TRUE), function(x) paste(x[-4], collapse = "\t"), "")
  })
  on.exit(unlink(v3))
  p = hmf_read_panel(v3)
  expect_equal(p[gene == "MYC", amp_ratio], 2)
  expect_equal(p[gene == "BRAF", amp_ratio], HMF_AMP_RATIO_DEFAULT)
  expect_false(any(p$report_splice))
})

test_that("known fusions: types, exon ranges and IG/ENHANCER overrides in both formats", {
  fu = hmf_read_fusions(file.path(hmf_fx, "known_fusion_data.38.csv"))
  expect_equal(attr(fu, "format"), "v2 (IG_*)")
  expect_setequal(fu$type, c("KNOWN_PAIR", "PROMISCUOUS_3", "EXON_DEL_DUP", "ENHANCER_KNOWN_PAIR",
                             "ENHANCER_PROMISCUOUS"))
  alk = fu[three_gene == "ALK"]
  expect_equal(c(alk$down_min, alk$down_max), c(3L, 4L))
  expect_true(alk$high_impact)
  expect_true(is.na(alk$five_gene))
  ig = fu[type == "ENHANCER_KNOWN_PAIR"]
  expect_equal(fu[type == "ENHANCER_PROMISCUOUS", .(five_gene, enh_orient, enh_chrom)],
               data.table(five_gene = "IGK", enh_orient = 0L, enh_chrom = "chr2"))
  expect_equal(list(ig$enh_orient, ig$enh_chrom, ig$enh_start, ig$enh_end, ig$down_distance),
               list(-1L, "chr14", 1000L, 5000L, 500L))

  v3 = hmf_v3_copy("known_fusion_data.38.csv",
                   function(l) gsub("IG_RANGE", "ENHANCER_RANGE", sub("^IG_", "ENHANCER_", l)))
  on.exit(unlink(v3))
  fu3 = hmf_read_fusions(v3)
  expect_equal(attr(fu3, "format"), "v3 (ENHANCER_*)")
  expect_equal(fu3[type == "ENHANCER_KNOWN_PAIR", enh_start], 1000L)
})

test_that("the Ensembl cache keeps canonical transcripts only, ranked, chr-prefixed", {
  ens = hmf_read_ensembl(file.path(hmf_fx, "ensembl_data"))
  expect_equal(nrow(ens$genes), 8L)
  expect_true(all(startsWith(ens$genes$chrom, "chr")))
  expect_false(any(startsWith(ens$exons$transcript, "ENST9")))
  expect_equal(ens$exons[gene == "TP53", rank], 1:4)
  expect_equal(ens$exons[gene == "TP53" & rank == 1, start], 20700L)   # minus strand
  expect_equal(nrow(hmf_read_ensembl(file.path(hmf_fx, "ensembl_data"), "MYC")$genes), 1L)
})

test_that("hotspots split INFO/input into gene and protein; the BED becomes 1-based", {
  hs = hmf_read_hotspots(file.path(hmf_fx, "KnownHotspots.somatic.38.vcf"))
  expect_equal(hs[pos == 10650, .(hs_gene, hs_protein)], data.table(hs_gene = "BRAF", hs_protein = "p.V600E"))
  expect_true(is.na(hs[hs_gene == "TERT", hs_protein]))
  bed = hmf_read_bed(file.path(hmf_fx, "ActionableCodingPanel.38.bed"))
  expect_equal(bed[gene == "BRAF", c(start, end)], c(10600L, 10700L))
})

# ---- Small variants --------------------------------------------------------------------

vep_rows = function() data.table(
  chrom = c("chr7", "chr17", "chr17", "chr2", "chr8", "chr5", "chr7"),
  pos   = c(10650L, 20350L, 20360L, 5350L, 50250L, 1295113L, 10661L),
  ref   = c("A", "C", "G", "T", "A", "G", "C"),
  alt   = c("T", "T", "A", "C", "G", "A", "A"),
  symbol = c("BRAF", "TP53", "TP53", "ALK", "MYC", NA, "BRAF"),
  gene_id = c("ENSG05", "ENSG06", "ENSG06", "ENSG03", "ENSG08", NA, "ENSG05"),
  consequence = c("missense_variant", "stop_gained", "missense_variant", "missense_variant",
                  "missense_variant", "intergenic_variant", "missense_variant"),
  hgvsp = c("ENSP05:p.Val600Glu", "ENSP06:p.Arg213Ter", "ENSP99:p.Gly1Ser", "ENSP03:p.Phe1174Leu",
            "ENSP08:p.Asn86Ser", NA, "ENSP05:p.Val601Lys"),
  # TP53's second row is annotated on a non-canonical transcript only
  feature = c("ENST05.1", "ENST06", "ENST906", "ENST03", "ENST08", NA, "ENST05"),
  canonical = c("YES", "YES", NA, "YES", "YES", NA, "YES"),
  coord_space = "vcf")

test_that("small variants follow the panel flags on the canonical transcript", {
  sm = hmf_small_variants(vep_rows(), NULL, hmf_res())
  braf = sm[gene == "BRAF" & pos == 10650]
  expect_true(braf$reported); expect_equal(braf$hotspot, "exact"); expect_equal(braf$reason, "hotspot")
  expect_true(braf$in_panel_region)
  tp53 = sm[gene == "TP53"]
  expect_equal(nrow(tp53), 1L)                       # the non-canonical row is not classified
  expect_equal(tp53$effect, "NONSENSE_OR_FRAMESHIFT"); expect_true(tp53$reported)
  expect_true(sm[gene == "ALK", reported])
  myc = sm[gene == "MYC"]
  expect_false(myc$reported); expect_equal(myc$reason, "reportMissense is false")
})

test_that("a hotspot outside the panel is listed once, by its own gene, not reported", {
  sm = hmf_small_variants(vep_rows(), NULL, hmf_res())
  tert = sm[pos == 1295113]
  expect_equal(nrow(tert), 1L)
  expect_equal(tert$gene, "TERT")
  expect_false(tert$reported); expect_false(tert$in_panel)
})

test_that("the same protein change as a hotspot with other alleles is a hotspot match", {
  # VEP's 3-letter p.Val601Lys against SERVE's 1-letter p.V601K
  sm = hmf_small_variants(vep_rows(), NULL, hmf_res())
  expect_equal(sm[pos == 10661, hotspot], "protein")
  expect_true(sm[pos == 10661, reported])
  expect_equal(sm[pos == 10661, reason], "hotspot (same protein change)")
  expect_equal(hmf_protein_1letter(c("p.Val600Glu", "p.Arg213Ter", "p.Gly12_Gly13delinsValCys", "p.V600E")),
               c("p.V600E", "p.R213*", "p.G12_G13delinsVC", "p.V600E"))
})

test_that("a protein hotspot match is not the reason when reportSomaticHotspot is false", {
  no_hs = hmf_v3_copy("DriverGenePanel.38.tsv", function(l) {
    i = startsWith(l, "BRAF\t"); x = strsplit(l[i], "\t", fixed = TRUE)[[1]]
    x[8] = "false"; l[i] = paste(x, collapse = "\t"); l
  })
  on.exit(unlink(no_hs))
  sm = hmf_small_variants(vep_rows(), NULL, hmf_res(driver_panel = no_hs))
  expect_equal(sm[pos == 10661, hotspot], "protein")
  expect_true(sm[pos == 10661, reported])                    # reportMissense
  expect_equal(sm[pos == 10661, reason], "missense reportable")
})

test_that("without a hotspot file, a hotspot-only call is no longer reported", {
  sm = hmf_small_variants(vep_rows(), NULL, hmf_res(hotspots = NULL))
  expect_equal(sm[gene == "BRAF" & pos == 10650, hotspot], "")
  expect_true(sm[gene == "BRAF" & pos == 10650, reported])   # still reportMissense
  expect_false(any(sm$pos == 1295113))
})

# ---- Copy number -----------------------------------------------------------------------

cn_seg = function(chrom, start, end, total, zero = FALSE) {
  data.table(chrom = chrom, start = start, end = end, cn_a = total, cn_b = 0, total = total,
             minor = 0, conf = NA_real_, sv_ids = NA_character_, conf_a = NA_real_,
             conf_b = NA_real_, zero_coverage = zero)
}

test_that("ASCAT and Wakhan are called separately; zero-coverage Wakhan segments are no data", {
  wk = rbind(cn_seg("chr8", 49000, 52000, 10), cn_seg("chr9", 29000, 32000, 0),
             cn_seg("chr17", 19000, 22000, 0, zero = TRUE))
  as = rbind(cn_seg("chr8", 49000, 52000, 5), cn_seg("chr9", 29000, 32000, 0),
             cn_seg("chr17", 19000, 22000, 0))
  cna = list(source = "wakhan", ploidy = 2, segments = wk,
             profiles = list(list(rank = 1L, ploidy = 2, segments = wk)),
             ascat_segments = as, ascat_ploidy = 2)
  cn = hmf_copy_number(cna, NULL, hmf_res())
  expect_equal(cn[gene == "MYC", source], "Wakhan")          # 10 > 3 x 2; ASCAT's 5 is not
  expect_equal(cn[gene == "MYC", call], "AMP")
  expect_setequal(cn[gene == "CDKN2A", source], c("Wakhan", "ASCAT"))
  expect_equal(unique(cn[gene == "CDKN2A", support]), "Wakhan, ASCAT")
  # TP53: deleted on ASCAT, uncovered (no reads) on Wakhan
  expect_equal(cn[gene == "TP53", source], "ASCAT")
  expect_true("TP53" %in% attr(cn, "uncovered")$wakhan)
  expect_equal(unname(attr(cn, "sources")), c("Wakhan (ploidy 2)", "ASCAT (ploidy 2)"))
})

test_that("an SV-caller deletion over reportDeletion exons is listed whatever the CN says", {
  sv = data.table(id = "d1", caller = "severus", svtype = "DEL", chrom_a = "chr9", pos_a = 30450L,
                  chrom_b = "chr9", pos_b = 30650L, vaf = 0.4, filter = "PASS")
  cn = hmf_copy_number(NULL, sv, hmf_res())
  expect_equal(nrow(cn), 1L)
  expect_equal(cn$source, "Severus SV")
  expect_match(cn$note, "exon 2 of 4")
  # BRAF has reportDeletion false: a deletion over it is not listed
  sv2 = copy(sv)[, `:=`(chrom_a = "chr7", chrom_b = "chr7", pos_a = 10000L, pos_b = 10800L)]
  expect_equal(nrow(hmf_copy_number(NULL, sv2, hmf_res())), 0L)
})

test_that("PARTIAL_AMP when only part of the gene exceeds the threshold; male chrX halves it", {
  seg = rbind(cn_seg("chr8", 49000, 50450, 2), cn_seg("chr8", 50451, 52000, 7))
  g = hmf_cn_calls(seg, 2, "ASCAT", hmf_res())
  expect_equal(g[gene == "MYC", call], "PARTIAL_AMP")
})

# ---- Fusions ---------------------------------------------------------------------------

sv_row = function(id, a_chrom, a_pos, a_or, b_chrom, b_pos, b_or, svtype = "BND")
  data.table(id = id, caller = "severus", svtype = svtype, chrom_a = a_chrom, pos_a = a_pos,
             orient_a = a_or, chrom_b = b_chrom, pos_b = b_pos, orient_b = b_or,
             gene_a = NA_character_, gene_b = NA_character_)

test_that("a known pair is placed on exons and phase-checked", {
  # EWSR1 (+) intron 2 kept upstream, ATF1 (+) intron 2 joined downstream -> exon 3
  fu = hmf_fusions(sv_row("b1", "chr22", 1350L, 1L, "chr12", 1350L, -1L), hmf_res())
  expect_equal(fu[, .(type, five_gene, three_gene, five_exon, three_exon, frame, likelihood, reported)],
               data.table(type = "KNOWN_PAIR", five_gene = "EWSR1", three_gene = "ATF1",
                          five_exon = "intron 2 → exon 2", three_exon = "intron 2 → exon 3",
                          frame = "inframe", likelihood = "HIGH", reported = TRUE))
  # The same SV written from the other side finds the same fusion
  fu2 = hmf_fusions(sv_row("b2", "chr12", 1350L, -1L, "chr22", 1350L, 1L), hmf_res())
  expect_equal(fu2$five_gene, "EWSR1")
  # Wrong orientation: ATF1 would be the 5' partner, which is no known pair
  expect_equal(nrow(hmf_fusions(sv_row("b3", "chr22", 1350L, -1L, "chr12", 1350L, 1L), hmf_res())), 0L)
})

test_that("a promiscuous 3' partner fuses with any gene, within its exon range", {
  # EML4 (+) intron 2 -> ALK (-) intron 2, so ALK's exon 3: inside the 3;4 range
  fu = hmf_fusions(sv_row("p1", "chr2", 1350L, 1L, "chr2", 5450L, 1L), hmf_res())
  expect_equal(fu$type, "PROMISCUOUS_3")
  expect_equal(fu$likelihood, "HIGH")
  expect_true(fu$high_impact)
  # ALK intron 1 -> exon 2, outside the range: LOW, and not high-impact
  fu2 = hmf_fusions(sv_row("p2", "chr2", 1350L, 1L, "chr2", 5650L, 1L), hmf_res())
  expect_equal(fu2$likelihood, "LOW"); expect_false(fu2$high_impact)
})

test_that("an exonic breakend cannot make a promiscuous fusion reportable", {
  # EML4 (+) exon 3 -> ALK (-) intron 2: no phase check on the exonic side
  fu = hmf_fusions(sv_row("p3", "chr2", 1450L, 1L, "chr2", 5450L, 1L), hmf_res())
  expect_equal(fu$type, "PROMISCUOUS_3")
  expect_equal(fu$frame, "exonic (unchecked)")
  expect_false(fu$reported)
})

test_that("a promiscuous enhancer needs a 3' gene beside its partner, and none inside the locus", {
  # IGK range chr2:2000-3000; the partner breakend 500 bp upstream of MYC (+)
  fu = hmf_fusions(sv_row("i1", "chr2", 2500L, 1L, "chr8", 49500L, -1L), hmf_res())
  expect_equal(fu[, .(type, five_gene, three_gene, reported)],
               data.table(type = "ENHANCER_PROMISCUOUS", five_gene = "IGK", three_gene = "MYC",
                          reported = FALSE))
  # No gene near the partner breakend, or both breakends inside the locus (V(D)J recombination)
  expect_equal(nrow(hmf_fusions(sv_row("i2", "chr2", 2500L, 1L, "chr8", 3000000L, -1L), hmf_res())), 0L)
  expect_equal(nrow(hmf_fusions(sv_row("i3", "chr2", 2200L, 1L, "chr2", 2800L, -1L, "DEL"), hmf_res())), 0L)
})

test_that("EXON_DEL_DUP needs both fused exons in the listed ranges", {
  # TP53 (-): a deletion from intron 2 to intron 1 removes exon 2 -> exon 1 joined to exon 3
  fu = hmf_fusions(sv_row("d1", "chr17", 20450L, 1L, "chr17", 20650L, -1L, "DEL"), hmf_res())
  expect_equal(fu$type, "EXON_DEL_DUP")
  expect_equal(c(fu$five_exon, fu$three_exon), c("intron 1 → exon 1", "intron 2 → exon 3"))
  # Removing exon 3 instead is not the listed event
  expect_equal(nrow(hmf_fusions(sv_row("d2", "chr17", 20250L, 1L, "chr17", 20450L, -1L, "DEL"), hmf_res())), 0L)
})

test_that("an IG enhancer known pair needs the IG breakend's orientation and the gene nearby", {
  fu = hmf_fusions(sv_row("e1", "chr14", 2000L, -1L, "chr8", 49000L, -1L), hmf_res())
  expect_equal(fu$type, "ENHANCER_KNOWN_PAIR")
  expect_equal(c(fu$five_gene, fu$three_gene), c("IGH", "MYC"))
  expect_true(fu$reported)
  expect_equal(nrow(hmf_fusions(sv_row("e2", "chr14", 2000L, 1L, "chr8", 49000L, -1L), hmf_res())), 0L)
  expect_equal(nrow(hmf_fusions(sv_row("e3", "chr14", 2000L, -1L, "chr8", 2000000L, -1L), hmf_res())), 0L)
})

test_that("no fusion file, or SVs without orientation, give an empty table", {
  sv = sv_row("b1", "chr22", 1350L, 1L, "chr12", 1350L, -1L)
  expect_equal(nrow(hmf_fusions(sv, hmf_res(known_fusions = NULL))), 0L)
  expect_equal(nrow(hmf_fusions(sv[, !c("orient_a", "orient_b")], hmf_res())), 0L)
})

# ---- SV orientation, entry point, gene list --------------------------------------------

test_that("breakend orientation follows the VCF bracket notation and the symbolic types", {
  o = .sv_orientation(c("N[chr2:5[", "N]chr2:5]", "]chr2:5]N", "[chr2:5[N", "N.", ".N",
                        "<DEL>", "<DUP>", "<INV>", "<INV>", "<INS>"),
                      c("BND", "BND", "BND", "BND", "sBND", "sBND", "DEL", "DUP", "INV", "INV", "INS"),
                      c(NA, NA, NA, NA, NA, NA, NA, NA, "++", "+-", NA))
  expect_equal(o$own,     c(1L, 1L, -1L, -1L, 1L, -1L, 1L, -1L, 1L, NA, NA))
  expect_equal(o$partner, c(-1L, 1L, 1L, -1L, NA, NA, -1L, 1L, 1L, NA, NA))
})

test_that("build_hmf is NULL when off and reports a load failure instead of erroring", {
  expect_null(build_hmf(NULL, NULL, NULL, list(), list()))
  bad = build_hmf(hmf_files(ensembl_data_dir = tempfile()), NULL, NULL, list(), list())
  expect_match(bad$error, "Ensembl data cache file not found")
  ok = build_hmf(hmf_files(), vep_rows(), NULL, list(), list(sex = "female"))
  expect_true(ok$small[gene == "BRAF" & pos == 10650, reported])
  expect_equal(ok$n_genes, 7L)
  expect_equal(ok$n_panel_columns, 15L)
})

test_that("the driver panel becomes a coordinate gene list with roles", {
  tsv = tempfile(fileext = ".tsv"); on.exit(unlink(tsv))
  hmf_gene_list_tsv(file.path(hmf_fx, "DriverGenePanel.38.tsv"), file.path(hmf_fx, "ensembl_data"), tsv)
  p = load_gene_panel(tsv, "hg38")
  expect_true(p$has_coords)
  expect_equal(p$reference, "hg38")
  expect_setequal(p$genes, c("BRAF", "TP53", "CDKN2A", "MYC", "ALK", "EWSR1", "ATF1"))
  roles = panel_roles(p)
  expect_equal(roles[match(c("TP53", "MYC"), gene), role], c("tsg", "oncogene"))
  expect_error(load_gene_panel(tsv, "t2t"), "declares reference")
})

# ---- CHM13: a lifted copy declares its build; VEP's transcript ids are another namespace ---

# A copy of the fixture Ensembl cache declaring `ref` in a REFERENCE file, as HMF_LIFTOVER writes it
hmf_lifted_ensembl = function(ref = "t2t") {
  d = tempfile("ensembl_lifted")
  dir.create(d)
  file.copy(list.files(file.path(hmf_fx, "ensembl_data"), full.names = TRUE), d)
  writeLines(ref, file.path(d, "REFERENCE"))
  d
}

test_that("on VEP ids outside HMF's namespace, VEP's canonical transcript stands in", {
  # The CHM13 rapid-release cache: ENSG0522.../ENST0522... ids, CANONICAL set, no MANE
  t2t = vep_rows()[, `:=`(gene_id = sub("^ENSG", "ENSG0522", gene_id),
                          feature = sub("^ENST", "ENST0522", feature))]
  sm = hmf_small_variants(t2t, NULL, hmf_res())
  expect_equal(attr(sm, "transcript_rule"), "vep_canonical")
  expect_true(sm[gene == "TP53", reported])                 # kept, matched by symbol
  expect_equal(nrow(sm[gene == "TP53"]), 1L)                # the non-canonical row still is not
  expect_true(sm[gene == "ALK", reported])
  # HMF's own ids: the HMF canonical transcript rule
  expect_equal(attr(hmf_small_variants(vep_rows(), NULL, hmf_res()), "transcript_rule"), "hmf_canonical")
})

test_that("the resources' build is declared by the Ensembl cache's REFERENCE file, else GRCh38", {
  expect_equal(hmf_declared_reference(hmf_files())$reference, "hg38")
  expect_null(hmf_check_reference(hmf_files(), "hg38"))
  expect_equal(hmf_check_reference(hmf_files(), "t2t")$kind, "undeclared")

  lifted = hmf_lifted_ensembl("T2T-CHM13v2.0")
  on.exit(unlink(lifted, recursive = TRUE))
  f = hmf_files(ensembl_data_dir = lifted)
  expect_equal(hmf_declared_reference(f)$reference, "t2t")
  expect_null(hmf_check_reference(f, "t2t"))
  expect_equal(hmf_check_reference(f, "hg38")$kind, "mismatch")

  # A hotspot VCF declaring another build than the cache: resources not lifted together
  hs = tempfile(fileext = ".vcf"); on.exit(unlink(hs), add = TRUE)
  l = readLines(file.path(hmf_fx, "KnownHotspots.somatic.38.vcf"))
  writeLines(c(l[1], "##assembly=GRCh38", l[-1]), hs)
  bad = hmf_check_reference(hmf_files(ensembl_data_dir = lifted, hotspots = hs), "t2t")
  expect_equal(bad$kind, "mismatch")
  expect_match(bad$msg, "lifted together")
})

test_that("the hmf_driver gene list carries the reference it was built for", {
  tsv = tempfile(fileext = ".tsv"); on.exit(unlink(tsv))
  hmf_gene_list_tsv(file.path(hmf_fx, "DriverGenePanel.38.tsv"), file.path(hmf_fx, "ensembl_data"), tsv, "t2t")
  expect_equal(load_gene_panel(tsv, "t2t")$reference, "t2t")
  expect_error(load_gene_panel(tsv, "hg38"), "declares reference")
})
